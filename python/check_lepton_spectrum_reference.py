#!/usr/bin/env python3
"""Independent single-muon bin actions, selected hard geometry and O(alpha).

Literal lattice pairs and direct allowed-azimuth integration replace the C++
CDF/prefix and spherical-cap actions. Selected hard matrices are sampled, not
fully reconstructed; inherited coupled tests audit the full contraction.
"""
import argparse
import bisect
import functools
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import source_digest
from check_xspace_reference import independent_cells, power_evolution, convolution, add, unit, scale
from check_coupled_decay_reference import hard_cells, causal
from check_fiducial_reference import gauss, band_quadrature
from check_moment_reference import analytic_running


def bounds(run, p):
    pt = run['pt_min']
    return min(math.tanh(run['eta_max']), math.sqrt((1-pt/p)*(1+pt/p))) if p > pt else 0.


def restrict(b, p, cut, threshold):
    if cut == 0:
        return b
    if cut == 2:
        return min(b, math.tanh(threshold))
    return min(b, math.sqrt((1-threshold/p)*(1+threshold/p))) if p > threshold else 0.


def angular(run, p, q, cosine, sine, cut, threshold):
    b, c = bounds(run, p), bounds(run, q)
    return (band_quadrature(restrict(b, p, cut, threshold), c, cosine, sine)
            +band_quadrature(b, restrict(c, q, cut, threshold), cosine, sine))/2


def born(run, f, g, cut, threshold, mass_cut=True):
    M = run['parent_mass']
    if mass_cut and not run['mass_min'] < M*math.sqrt(f*g) < run['mass_max']:
        return 0.
    p, q = M*f/2, M*g/2
    base = min(bounds(run, p), bounds(run, q))
    return (restrict(base, p, cut, threshold)+restrict(base, q, cut, threshold))/2


def panels(function, lo, hi, edges, count):
    if lo >= hi:
        return 0.
    edges = sorted(set([lo, hi]+[v for v in edges if lo < v < hi]))
    value = 0.
    for a, b in zip(edges, edges[1:]):
        def mapped(u):
            s, c = math.sin(u), math.cos(u)
            return function(a+(b-a)*s*s)*2*(b-a)*s*c
        value += gauss(mapped, 0, math.pi/2, count)
    return value


def energies(run, cut, t):
    output = [run['pt_min'], run['pt_min']*math.cosh(run['eta_max'])]
    if cut == 1:
        output.extend((t, t*math.cosh(run['eta_max'])))
    elif cut == 2:
        output.append(run['pt_min']*math.cosh(t))
    return output


def angular_panels(function, lower, upper, breaks, transitions, count):
    # Besides energy thresholds, split where the two spherical caps become
    # nested, separate, or touch across the opposite direction. The allowed-
    # azimuth integrand has fractional-power contacts at these boundaries.
    edges = sorted(set([lower, upper]+[v for v in breaks if lower < v < upper]))
    contacts = list(edges)
    for lo, hi in zip(edges, edges[1:]):
        steps = max(8, count//6)
        probes = [lo+(hi-lo)*k/steps for k in range(steps+1)]
        values = [transitions(u) for u in probes]
        for left, right, va, vb in zip(probes, probes[1:], values, values[1:]):
            for k, (fa, fb) in enumerate(zip(va, vb)):
                if fa == 0:
                    contacts.append(left)
                if fb == 0:
                    contacts.append(right)
                if fa*fb >= 0:
                    continue
                a, b = left, right
                for _ in range(48):
                    mid = (a+b)/2
                    fm = transitions(mid)[k]
                    if fa*fm <= 0:
                        b = mid
                    else:
                        a, fa = mid, fm
                contacts.append((a+b)/2)
    resolved = [lower]
    for u in sorted(contacts):
        if u-resolved[-1] > 1e-12 and upper-u > 1e-12:
            resolved.append(u)
    return panels(function, lower, upper, resolved[1:], max(16, count//2))


def geometry(run, row, f, g, cut, threshold, photon, count):
    M, pt = run['parent_mass'], run['pt_min']
    if min(f, g) <= 2*pt/M:
        return 0.
    lo = (run['mass_min']/M)**2/(f*g)
    hi = min(1. if photon else 1-row['soft_cutoff'], (run['mass_max']/M)**2/(f*g))
    E = energies(run, cut, threshold)
    outer_edges = [2*e/(M*z) for e in E for z in (f, g)]
    outer_edges += [4*e/(M*z)-1 for e in E for z in (f, g)]
    outer_edges.extend((2*pt/M*(1/f+1/g)-1, g/f))
    if not photon:
        outer_edges.append(f/g)  # second collinear orientation changes branch here
    if cut == 1:
        outer_edges.extend((2*(pt+threshold)/M-1, g/f*pt/threshold, g/f*threshold/pt))
        outer_edges.extend((2*(threshold/f+pt/g)/M-1, 2*(pt/f+threshold/g)/M-1))
    umin = math.asin(math.sqrt(row['collinear_cutoff']))

    def outer(z):
        reference = born(run, z*f, g, cut, threshold)
        reverse = born(run, z*g, f, cut, threshold)
        inner_edges = []
        if photon:
            contacts = [2*e/(M*f)-z for e in E]+[1-2*e/(M*g) for e in E]
            inner_edges = [math.asin(math.sqrt(t/(1-z))) for t in contacts if 0 < t < 1-z]
        else:
            for e in E:
                for fraction in (f, g):
                    c = abs(((1+z)-4*e/(M*fraction))/(1-z))
                    if 0 < c < 1:
                        inner_edges.append(math.asin(math.sqrt(1-c)))

        def integrand(u):
            s, c = math.sin(u), math.cos(u)
            if photon:
                t = (1-z)*s*s
                x, y = z+t, 1-t
                cosine = min(1., abs(1-2*z/(x*y)))
                sine = 2*math.sqrt(max(0., z*t*(1-z-t)))/(x*y)
                accepted = angular(run, M*x*f/2, M*y*g/2, cosine, sine, cut, threshold)
                density = ((1-z)**2+y*y)/(z*t)
                return density*(accepted-reference)*2*(1-z)*s*c
            cosine_angle = 1-s*s
            x, y = (1+z-(1-z)*cosine_angle)/2, (1+z+(1-z)*cosine_angle)/2
            cosine = min(1., abs(x*y-2*z)/(x*y))
            sine = math.sqrt(z)*(1-z)*math.sqrt((1-cosine_angle)*(1+cosine_angle))/(x*y)
            a = angular(run, M*x*f/2, M*y*g/2, cosine, sine, cut, threshold)
            b = angular(run, M*x*g/2, M*y*f/2, cosine, sine, cut, threshold)
            density = (x*x+y*y)/((1-x)*(1-y))*(1-z)/2
            return density*(a-reference+b-reverse)*2*s*c

        def transitions(u):
            s, c = math.sin(u), math.cos(u)
            if photon:
                t = (1-z)*s*s
                x, y = z+t, 1-t
                cosine = min(1., abs(1-2*z/(x*y)))
                sine = 2*math.sqrt(max(0., z*t*(1-z-t)))/(x*y)
            else:
                a = 1-s*s
                x, y = (1+z-(1-z)*a)/2, (1+z+(1-z)*a)/2
                cosine = min(1., abs(x*y-2*z)/(x*y))
                sine = math.sqrt(z)*(1-z)*math.sqrt((1-a)*(1+a))/(x*y)
            angle = math.atan2(sine, cosine)
            result = []
            fractions = [(f, g)] if photon or f == g else [(f, g), (g, f)]
            for first, second in fractions:
                p, q = M*x*first/2, M*y*second/2
                b, c = bounds(run, p), bounds(run, q)
                pairs = [(b, c)] if cut == 0 else [(restrict(b, p, cut, threshold), c),
                                                   (b, restrict(c, q, cut, threshold))]
                for one, two in pairs:
                    alpha, beta = math.acos(one), math.acos(two)
                    result.extend((angle-abs(alpha-beta), angle-alpha-beta, math.pi-angle-alpha-beta))
            return result

        return angular_panels(integrand, umin, math.pi/2, inner_edges, transitions, count)

    return panels(outer, lo, hi, outer_edges, count)


def projected(run, cut, t):
    M, pt = run['parent_mass'], run['pt_min']
    B = born(run, 1., 1., cut, t)
    L = 2*math.log(M/run['muon_mass'])
    edges = [2*pt/M, (run['mass_min']/M)**2, (run['mass_max']/M)**2]
    edges += [2*e/M for e in energies(run, cut, t)]
    if cut == 1:
        edges.append(pt/t)
    def integrand(z):
        return 2*(1+z*z)/(1-z)*(L+math.log(z)-1)*(born(run, z, 1., cut, t)-B)
    return 1.5*B+panels(integrand, 0., 1., edges, 96)


def bin_values(run, observable, function):
    edges = run[observable+'_edges']
    if observable == 'pt':
        def cumulative(t):
            return 0. if t >= run['parent_mass']/2 else function(1, t)
        return [2*(cumulative(a)-cumulative(b)) for a, b in zip(edges, edges[1:])]
    def cumulative(t):
        return 0. if t == 0 else function(2, t)
    return [(cumulative(-a)-cumulative(-b) if b <= 0 else
             cumulative(b)-cumulative(a) if a >= 0 else cumulative(-a)+cumulative(b))
            for a, b in zip(edges, edges[1:])]


def verify(run):
    row = run['points'][0]
    if run['mu_0'] != run['muon_mass'] or run['mu_frag'] != run['parent_mass']:
        raise ValueError('Independent source fixture requires mu_0=m_mu and mu_frag=M')
    n, span = row['bins'], row['log_span']
    h, a, length = span/n, row['a_initial'], len(row['same'])
    source = independent_cells(n, span, 64)
    p0, p1, b1 = [[[[v for v in channel[:length]] for channel in parent] for parent in A] for A in source]
    L = 2*math.log(run['mu_frag']/run['mu_0'])
    D = convolution(power_evolution(p0, p1, a, L, 20), add(unit(length), b1, a))
    S, O, G = D[0][0], D[0][1], D[2][0]
    c1 = add(b1, p0, L)
    c2 = add(scale(p1, L), convolution(p0, b1), L)
    c2 = add(c2, convolution(p0, p0), L*L/2)
    c2 = add(c2, p0, L*L/3)
    C, CO, CG, C2 = c1[0][0], c1[0][1], c1[2][0], c2[0][0]
    I = [1.]+[0.]*(length-1)
    hm, hg = [v[:length] for v in hard_cells(n, span, 64)]
    errors = {'ff_weights': max(abs(x-y) for key, values in (('same', S), ('opposite', O), ('gamma', G))
                               for x, y in zip(row[key], values)),
              'hard_projection': max(abs(x-y) for key, values in (('hard_muon', hm), ('hard_photon', hg))
                                     for x, y in zip(row[key], values)),
              'binned_ff': 0., 'binned_muon_projection': 0., 'binned_photon_projection': 0.,
              'selected_geometry': 0., 'geometry_refinement': 0., 'continuum_coefficient1': 0.,
              'continuum_refinement': 0., 'fixed_order': 0., 'born_bins': 0.,
              'endpoint_geometry': 0., 'endpoint_refinement': 0., 'compression_moments': 0.,
              'forbidden_support': 0., 'binned_muon_geometry_algebra': 0., 'binned_photon_geometry_algebra': 0.,
              'lattice_coefficient1': 0., 'coefficient2_algebra': 0., 'coefficient2_azimuth': 0.,
              'coefficient2_azimuth_refinement': 0., 'matched_algebra': 0.}
    plus = lambda x, y: [u+v for u, v in zip(x, y)]

    def literal(u, v, observable):
        output = [0.]*(len(run[observable+'_edges'])-1)
        M = run['parent_mass']
        for i, x in enumerate(u):
            for j, y in enumerate(v):
                f, g = math.exp(-h*i), math.exp(-h*j)
                if i+j == 0:
                    fraction = float(run['mass_min'] < M < run['mass_max'])
                else:
                    center = (i+j)*h
                    loss_lo, loss_hi = 2*math.log(M/run['mass_max']), 2*math.log(M/run['mass_min'])
                    fraction = max(0., min(center+h/2, loss_hi)-max(center-h/2, loss_lo))/h
                angular_bins = bin_values(run, observable, lambda cut, t: born(run, f, g, cut, t, False))
                for k, value in enumerate(angular_bins):
                    output[k] += x*y*value*fraction/(f*g)
        return output

    @functools.lru_cache(None)
    def selected_geometry(f, g, cut, t, photon, count):
        return geometry(run, row, f, g, cut, t, photon, count)

    operators = {}
    shared = None
    minimum = max(2*run['pt_min']/run['parent_mass'],
                  (run['mass_min']/run['parent_mass'])**2*math.exp(-h/2))
    loss = math.sqrt(1-minimum)
    contacts = [loss*k/row['ff_panels'] for k in range(1, row['ff_panels'])]
    contacts += [math.sqrt(1-z) for z in run['ff_fraction_breaks'] if z > minimum]
    loss_edges = [0.]
    margin = 1024*2.**(1-run['real_digits'])
    for value in sorted(contacts):
        if value-loss_edges[-1] > margin and loss-value > margin:
            loss_edges.append(value)
    loss_edges.append(loss)
    def compress(weights, op):
        f = op['fractions']
        output = [weights[0]]+[0.]*(len(f)-1)
        for j, weight in enumerate(weights[1:], 1):
            z = math.exp(-h*j)
            k = min(len(loss_edges)-2, bisect.bisect_right(loss_edges, math.sqrt(1-z))-1)
            hi, lo = op['panel_edges'][k:k+2]
            indices = [i for i, x in enumerate(f) if i and lo < x < hi]
            if len(indices) != 4:
                raise RuntimeError('Nonuniform panel lost its cubic nodes')
            for i in indices:
                output[i] += weight/z*math.prod((z-f[k])/(f[i]-f[k]) for k in indices if k != i)
        for power in range(4):
            errors['compression_moments'] = max(errors['compression_moments'], abs(
                sum(w*f[i]**power for i, w in enumerate(output))
                -sum(w*math.exp(-h*j*(power-1)) for j, w in enumerate(weights))))
        return output

    def mixed(A, B, qa, qb, matrix, left, right):
        # Separate literal endpoint strips from the strictly interior tensor.
        return (A[0]*B[0]*matrix[0][0]
                +sum(A[j]*B[0]*left[j]+A[0]*B[j]*right[j] for j in range(1, length))
                +sum(qa[i]*matrix[i][j]*qb[j] for i in range(1, len(qa)) for j in range(1, len(qb))))

    for observable in ('pt', 'eta'):
        for op in row[observable+'_cumulatives']:
            key = op['cut'], op['threshold']
            if shared is None:
                shared = op['fractions'], op['panel_edges']
            elif shared != (op['fractions'], op['panel_edges']):
                raise RuntimeError('Cumulative measurements changed the common FF partition')
            if key in operators:
                continue
            qS, qO, qG = (compress(weights, op) for weights in (S, O, G))
            muon, first, second = (op[k] for k in ('endpoint_muon', 'endpoint_photon_first', 'endpoint_photon_second'))
            gm, gg = op['geometry_muon'], op['geometry_photon']
            direct = {}
            for count in (48, 96):
                values = []
                for photon, reverse, field in ((False, False, muon), (True, False, first), (True, True, second)):
                    current = []
                    for j in range(length):
                        z = math.exp(-h*j)
                        value = selected_geometry(1. if reverse else z, z if reverse else 1., *key, photon, count)
                        current.append(value/z)
                        if count == 96:
                            errors['endpoint_geometry'] = max(errors['endpoint_geometry'], abs(value-z*field[j]))
                    values.append(current)
                direct[count] = values
            errors['endpoint_refinement'] = max(errors['endpoint_refinement'],
                max(abs(u-v)*math.exp(-h*j) for fine, coarse in zip(direct[96], direct[48])
                    for j, (u, v) in enumerate(zip(fine, coarse))))
            if key[0] == 1:
                critical = 2*key[1]/run['parent_mass']
                for j in range(1, length):
                    if math.exp(-h*j) <= critical:
                        basis = compress([float(k == j) for k in range(length)], op)
                        errors['forbidden_support'] = max(errors['forbidden_support'],
                            max((abs(w) for f, w in zip(op['fractions'], basis) if f > critical), default=0.))
            operators[key] = {
                'muon_geometry': mixed(S, S, qS, qS, gm, muon, muon)+mixed(O, O, qO, qO, gm, muon, muon),
                'photon_geometry': 2*mixed(G, plus(S, O), qG, plus(qS, qO), gg, first, second),
                'gm00': gm[0][0], 'muon_C': sum(x*y for x, y in zip(C, muon)),
                'photon_C': sum(x*y for x, y in zip(CG, first)),
                **{f'gm00_{n}': direct[n][0][0] for n in (48, 96)},
                **{f'muon_C_{n}': sum(x*y for x, y in zip(C, direct[n][0])) for n in (48, 96)},
                **{f'photon_C_{n}': sum(x*y for x, y in zip(CG, direct[n][1])) for n in (48, 96)}}

    def action(cut, t, field):
        key = (0, 0) if (cut == 1 and t <= run['pt_min']) or (cut == 2 and t >= run['eta_max']) else (cut, t)
        return operators[key][field]

    for observable in ('pt', 'eta'):
        for op in row[observable+'_cumulatives']:
            for sample in op['geometry_samples']:
                f, g = (op['fractions'][sample[key]] for key in ('i', 'j'))
                for photon, field in ((False, 'muon'), (True, 'photon')):
                    coarse = selected_geometry(f, g, op['cut'], op['threshold'], photon, 48)
                    fine = selected_geometry(f, g, op['cut'], op['threshold'], photon, 96)
                    errors['selected_geometry'] = max(errors['selected_geometry'], abs(fine-sample[field]))
                    errors['geometry_refinement'] = max(errors['geometry_refinement'], abs(fine-coarse))
        ff = plus(literal(S, S, observable), literal(O, O, observable))
        muon = plus(literal(causal(hm, S), S, observable), literal(causal(hm, O), O, observable))
        photon = [2*x for x in literal(causal(hg, G), plus(S, O), observable)]
        geometric = {field: bin_values(run, observable, lambda cut, t: action(cut, t, field))
                     for field in ('muon_geometry', 'photon_geometry', 'gm00', 'muon_C', 'photon_C',
                                   'gm00_48', 'gm00_96', 'muon_C_48', 'muon_C_96', 'photon_C_48', 'photon_C_96')}
        first_born = literal(C, I, observable)
        hard_born = literal(hm, I, observable)
        second_born = [2*x+y+z for x, y, z in zip(literal(C2, I, observable),
                                                literal(C, C, observable), literal(CO, CO, observable))]
        second_muon = plus(literal(causal(hm, C), I, observable), literal(hm, C, observable))
        second_photon = [2*x for x in literal(causal(hg, CG), I, observable)]
        @functools.lru_cache(None)
        def coefficient(cut, t, count):
            return projected(run, cut, t)+selected_geometry(1., 1., cut, t, False, count)
        fine = bin_values(run, observable, lambda cut, t: coefficient(cut, t, 96))
        coarse = bin_values(run, observable, lambda cut, t: coefficient(cut, t, 48))
        born_bins = bin_values(run, observable, lambda cut, t: born(run, 1., 1., cut, t))
        for k, b in enumerate(row[observable]['spectrum']):
            result = b['integral']
            for key, values, target in (('binned_ff', ff, 'ff'), ('binned_muon_projection', muon, 'muon_projected'),
                                        ('binned_photon_projection', photon, 'photon_projected')):
                errors[key] = max(errors[key], abs(values[k]-result[target]))
            for field in ('muon_geometry', 'photon_geometry'):
                name = 'binned_'+field+'_algebra'
                errors[name] = max(errors[name], abs(geometric[field][k]-result[field]))
            first = 2*first_born[k]+hard_born[k]+geometric['gm00'][k]
            errors['lattice_coefficient1'] = max(errors['lattice_coefficient1'], abs(first-result['lattice_coefficient1']))
            second_base = second_born[k]+second_muon[k]+second_photon[k]+(2/3)*L*hard_born[k]
            def second_coefficient(suffix):
                return (second_base+2*geometric['muon_C'+suffix][k]+2*geometric['photon_C'+suffix][k]
                        +(2/3)*L*geometric['gm00'+suffix][k])
            errors['coefficient2_algebra'] = max(errors['coefficient2_algebra'], abs(second_coefficient('')-result['coefficient2']))
            errors['coefficient2_azimuth'] = max(errors['coefficient2_azimuth'], abs(second_coefficient('_96')-result['coefficient2']))
            errors['coefficient2_azimuth_refinement'] = max(errors['coefficient2_azimuth_refinement'],
                abs(second_coefficient('_96')-second_coefficient('_48')))
            matched = ff[k]+row['a_hard']*(muon[k]+photon[k]+geometric['muon_geometry'][k]+geometric['photon_geometry'][k])
            matched += a*(result['continuum_coefficient1']-first)
            errors['matched_algebra'] = max(errors['matched_algebra'], abs(matched-result['matched']))
            errors['born_bins'] = max(errors['born_bins'], abs(born_bins[k]-result['born']))
            errors['continuum_coefficient1'] = max(errors['continuum_coefficient1'], abs(fine[k]-result['continuum_coefficient1']))
            errors['continuum_refinement'] = max(errors['continuum_refinement'], abs(fine[k]-coarse[k]))
            errors['fixed_order'] = max(errors['fixed_order'], abs(born_bins[k]+a*fine[k]-result['fixed_order']))
    limits = {'ff_weights': 2e-9, 'hard_projection': 2e-7, 'binned_ff': 2e-9,
              'binned_muon_projection': 2e-7, 'binned_photon_projection': 2e-9,
              'selected_geometry': 1e-4, 'geometry_refinement': 2e-4,
              'continuum_coefficient1': 1e-5, 'continuum_refinement': 2e-5,
              'fixed_order': 2e-8, 'born_bins': 1e-13,
              'endpoint_geometry': 1e-4, 'endpoint_refinement': 2e-4, 'compression_moments': 2e-12,
              'forbidden_support': 1e-14, 'binned_muon_geometry_algebra': 2e-8, 'binned_photon_geometry_algebra': 2e-9,
              'lattice_coefficient1': 2e-7, 'coefficient2_algebra': 2e-5, 'coefficient2_azimuth': 1e-2,
              'coefficient2_azimuth_refinement': 3e-2, 'matched_algebra': 1e-10}
    print(json.dumps(errors, indent=2), flush=True)
    for key, value in errors.items():
        if not math.isfinite(value) or value > limits[key]:
            raise RuntimeError(f'Independent single-muon check failed: {key}')
    return errors


def scan_checks(run):
    points = {p['label']: p for p in run['points']}
    for point in points.values():
        for observable in ('pt', 'eta'):
            for key, value in point['total'].items():
                if abs(sum(b['integral'][key] for b in point[observable]['spectrum'])-2*value) > 1e-10:
                    raise RuntimeError(f'Two-entry closure failed: {point["label"]} {observable} {key}')
            for b in point[observable]['spectrum']:
                density = b['matched_density_per_GeV' if observable == 'pt' else 'matched_density_per_eta']
                if abs(density*(b['upper']-b['lower'])-b['integral']['matched']) > 1e-14:
                    raise RuntimeError('Density normalisation failed')
        if max(point['weighted_wrap_bound'], point['weighted_tail_bound']) > 1e-10:
            raise RuntimeError('FF domain bounds failed')
        if point['compression_moment_residual'] > 1e-12:
            raise RuntimeError('FF compression moment check failed')
    if 'hard_slivers' not in points:
        return {}, {}
    for label in ('ff_panels_x2', 'grid_x2', 'grid_x4', 'span_x2', 'ff_panels_x4'):
        if label not in points:
            continue
        for observable in ('pt', 'eta'):
            if any(abs(b['integral']['fixed_order']-c['integral']['fixed_order']) > 1e-12
                   for b, c in zip(points[label][observable]['spectrum'], points['baseline'][observable]['spectrum'])):
                raise RuntimeError('Grid/panel refinement changed continuum fixed order')
            if label.startswith('ff_panels') and any(
                    abs(b['integral']['coefficient2']-c['integral']['coefficient2']) > 2e-8
                    for b, c in zip(points[label][observable]['spectrum'], points['baseline'][observable]['spectrum'])):
                raise RuntimeError('Endpoint-corrected C2 retained panel dependence')
    shifts, variations = {}, {}
    comparisons = [('ff_panels', 'ff_panels_x2', 'baseline'),
                               ('grid_x2', 'grid_x2', 'ff_panels_x2'), ('grid_x4', 'grid_x4', 'grid_x2'),
                               ('span', 'span_x2', 'grid_x4'), ('hard_slivers', 'hard_slivers', 'span_x2'),
                   *((v, v, 'hard_slivers') for v in ('scale_half', 'scale_double', 'start_half', 'start_double'))]
    if 'ff_panels_x4' in points:
        comparisons.append(('ff_panels_further', 'ff_panels_x4', 'ff_panels_x2'))
    for label, left, right in comparisons:
        output = {'pair_total_shift': points[left]['total']['matched']-points[right]['total']['matched']}
        for observable in ('pt', 'eta'):
            delta = [b['integral']['matched']-c['integral']['matched']
                     for b, c in zip(points[left][observable]['spectrum'], points[right][observable]['spectrum'])]
            output[observable] = {'bins': delta, 'maximum_bin_shift': max(abs(x) for x in delta)}
        (variations if label.startswith(('scale_', 'start_')) else shifts)[label] = output
    for label in variations:
        point, central = points[label], points['hard_slivers']
        if abs(point['a_hard']-central['a_hard']) > 1e-13:
            raise RuntimeError('Scale variation changed hard coupling')
        for observable in ('pt', 'eta'):
            if any(abs(b['integral']['continuum_coefficient1']-c['integral']['continuum_coefficient1']) > 5e-6
                   for b, c in zip(point[observable]['spectrum'], central[observable]['spectrum'])):
                raise RuntimeError('Scale variation changed continuum first coefficient')
        if label.startswith('start_'):
            expected = analytic_running(run['alpha_msbar_at_mu0']/(2*math.pi),
                                        2*math.log(point['mu_0']/run['mu_0']), 'two_loop')[0]
            if abs(expected-point['a_initial']) > 1e-14:
                raise RuntimeError('Starting-scale coupling transport failed')
    return shifts, variations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('config', type=Path)
    parser.add_argument('bin_config', type=Path)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--panel-input', type=Path, help='Additional doubled-panel run at the baseline grid')
    parser.add_argument('--checked-coarse', type=Path, help='Reuse an audit with matching source/executable/input hashes')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.panel_input and not args.input:
        parser.error('--panel-input requires a completed --input scan')
    exe, config, bins = (v.resolve() for v in (args.executable, args.config, args.bin_config))
    root = Path(__file__).resolve().parents[1]
    if args.checked_coarse:
        coarse = json.loads(args.checked_coarse.read_text())
        hashes = coarse['provenance']
        if hashes['source_sha256'] != source_digest(root):
            raise RuntimeError('Checked coarse source is stale')
        for key, path in (('executable_sha256', exe), ('config_sha256', config), ('bin_config_sha256', bins)):
            if hashes[key] != hashlib.sha256(path.read_bytes()).hexdigest():
                raise RuntimeError(f'Checked coarse hash differs: {key}')
        if len(coarse['points']) != 1 or coarse['points'][0]['bins'] != 256 or coarse['points'][0]['ff_panels'] != 1:
            raise RuntimeError('Checked coarse artifact is not the standalone audit fixture')
        errors = coarse['independent_checks']
        scan_checks(coarse)
    else:
        coarse = json.loads(subprocess.check_output([str(exe), str(config), str(bins), '256', '1', '--dump'], text=True))
        errors = verify(coarse)
    run = json.loads(args.input.read_text()) if args.input else coarse
    input_keys = ('version', 'compiler', 'parent_mass', 'muon_mass', 'mu_0', 'mu_frag', 'alpha_msbar_at_mu0',
                  'pt_min', 'eta_max', 'mass_min', 'mass_max', 'pt_edges', 'eta_edges', 'ff_fraction_breaks', 'geometry_quadrature')
    for key in input_keys:
        if run[key] != coarse[key]:
            raise RuntimeError(f'Production/coarse inputs differ: {key}')
    if args.panel_input:
        extra = json.loads(args.panel_input.read_text())
        if any(extra[key] != run[key] for key in input_keys):
            raise RuntimeError('Further-panel inputs differ')
        reference = next(p for p in run['points'] if p['label'] == 'ff_panels_x2')
        point = extra['points'][0]
        for key in ('bins', 'log_span', 'soft_cutoff', 'collinear_cutoff', 'geometry_abs_tol',
                    'geometry_rel_tol', 'mu_frag', 'mu_0', 'alpha_msbar_at_mu0'):
            if point[key] != reference[key]:
                raise RuntimeError(f'Further-panel controls differ: {key}')
        if point['ff_panels'] != 2*reference['ff_panels']:
            raise RuntimeError('Further-panel run must double ff_panels_x2')
        point['label'] = 'ff_panels_x4'
        run['points'].append(point)
    run['independent_checks'] = errors
    run['refinement_rate_shifts'], run['variation_rate_shifts'] = scan_checks(run)
    if args.input:
        run['independent_coarse'] = coarse
    run['provenance'] = {'source_sha256': source_digest(Path(__file__).resolve().parents[1]),
        'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
        'config_sha256': hashlib.sha256(config.read_bytes()).hexdigest(),
        'bin_config_sha256': hashlib.sha256(bins.read_bytes()).hexdigest(),
        'input_json_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest() if args.input else None,
        'panel_input_json_sha256': hashlib.sha256(args.panel_input.read_bytes()).hexdigest() if args.panel_input else None,
        'checked_coarse_json_sha256': hashlib.sha256(args.checked_coarse.read_bytes()).hexdigest() if args.checked_coarse else None,
        'python_version': platform.python_version(), 'platform': platform.platform()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(run, indent=2, allow_nan=False)+'\n')
    print('Independent single-muon bin, selected hard geometry, Oalpha and two-entry closure checks passed')


if __name__ == '__main__':
    main()
