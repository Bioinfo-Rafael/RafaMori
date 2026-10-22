"""Balanced separated-reactant electronic and RRHO thermochemical differences."""
import numpy as np

KCAL = 627.509474


def validate_reference(reactant, ts, product, lys):
    """Reject inconsistent reference energies before calculating a barrier."""
    records = (reactant, ts, product, lys)
    for key in ('method_basis', 'solvent', 'multiplicity'):
        if any(r[key] != ts[key] for r in records):
            raise ValueError(f'Calculation mismatch: {key}')
    if 'UNKNOWN' in ts['method_basis'] or ts['solvent'] == 'UNKNOWN':
        raise ValueError('Unknown calculation conditions')
    if any(not r['normal_termination'] or not r['optimization_completed'] for r in records):
        raise ValueError('Optimization or normal termination missing')
    if reactant['charge'] + lys['charge'] != ts['charge'] or product['charge'] != ts['charge']:
        raise ValueError('Unbalanced charge')
    composition = reactant['composition'].copy()
    for z, count in lys['composition'].items():
        composition[z] = composition.get(z, 0) + count
    if composition != ts['composition'] or composition != product['composition']:
        raise ValueError('Unbalanced composition')
    for key in ('temperature_K', 'pressure_atm'):
        if any(r[key] is None or not np.isclose(r[key], ts[key]) for r in records):
            raise ValueError(f'Thermochemistry mismatch: {key}')
    for r in records:
        if r['G_hartree'] is None or r['frequency_count'] != 3*r['atom_count']-6:
            raise ValueError('Complete nonlinear-molecule frequencies/thermochemistry required')
    if lys['imaginary_frequencies']:
        raise ValueError('Lys reference is not a verified minimum')


def energy_terms(reactant, ts, product, lys, ealpha, eleft):
    validate_reference(reactant, ts, product, lys)
    er = reactant['energy_hartree'] + lys['energy_hartree']
    gr = reactant['G_hartree'] + lys['G_hartree']
    values = dict(
        TS_E_hartree=ts['energy_hartree'], substrate_E_hartree=reactant['energy_hartree'],
        product_E_hartree=product['energy_hartree'], Lys_E_hartree=lys['energy_hartree'],
        TS_G_hartree=ts['G_hartree'], substrate_G_hartree=reactant['G_hartree'],
        product_G_hartree=product['G_hartree'], Lys_G_hartree=lys['G_hartree'],
        separated_reactants_E_hartree=er, separated_reactants_G_hartree=gr,
        absolute_electronic_barrier_kcal_mol=(ts['energy_hartree']-er)*KCAL,
        absolute_RRHO_G_barrier_kcal_mol=(ts['G_hartree']-gr)*KCAL,
        electronic_reaction_energy_kcal_mol=(product['energy_hartree']-er)*KCAL,
        RRHO_G_reaction_energy_kcal_mol=(product['G_hartree']-gr)*KCAL,
        TS_minus_substrate_E_hartree=ts['energy_hartree']-reactant['energy_hartree'],
        product_minus_substrate_E_hartree=product['energy_hartree']-reactant['energy_hartree'],
        TS_minus_substrate_G_hartree=ts['G_hartree']-reactant['G_hartree'],
        absolute_barrier_valid=True, reference_conditions_match=True,
        reactant_minimum_verified=not reactant['imaginary_frequencies'],
        RRHO_G_provisional=bool(reactant['imaginary_frequencies']),
        temperature_K=ts['temperature_K'], pressure_atm=ts['pressure_atm'],
        IRC_left_to_TS_kcal_mol=(ts['energy_hartree']-eleft)*KCAL,
        IRC_left_to_alpha_kcal_mol=(ealpha-eleft)*KCAL,
        alpha_to_TS_kcal_mol=(ts['energy_hartree']-ealpha)*KCAL,
        separated_reactants_to_alpha_kcal_mol=(ealpha-er)*KCAL,
        alpha_minus_substrate_E_hartree=ealpha-reactant['energy_hartree'])
    for key, label in [('H_hartree','enthalpy'), ('ZPE_corrected_hartree','ZPE_corrected_energy')]:
        values[f'{label}_barrier_kcal_mol']=(ts[key]-reactant[key]-lys[key])*KCAL
    values.update(entropy_terms(reactant,ts,product,lys))
    values.update(reverse_barrier_terms(ts, product))
    return values


def reverse_barrier_terms(ts, product):
    """Product-to-TS differences; independent of the separated reactants."""
    values = {
        f'reverse_{label}_barrier_kcal_mol': (ts[key] - product[key]) * KCAL
        for key, label in [('energy_hartree', 'electronic'),
                           ('H_hartree', 'enthalpy'), ('G_hartree', 'RRHO_G')]
    }
    tds = values['reverse_enthalpy_barrier_kcal_mol'] - values['reverse_RRHO_G_barrier_kcal_mol']
    values['T_delta_S_reverse_activation_kcal_mol'] = tds
    values['minus_T_delta_S_reverse_activation_kcal_mol'] = -tds
    values['reverse_RRHO_G_provisional'] = bool(product['imaginary_frequencies'])
    return values


def entropy_terms(reactant, ts, product, lys):
    """RRHO entropy and thermal corrections; electronic E is not enthalpy H."""
    records={'substrate':reactant,'TS':ts,'product':product,'Lys':lys}
    for record in records.values():
        if any(record.get(k) is None for k in ('energy_hartree','H_hartree','G_hartree','temperature_K')):
            raise ValueError('E, H, G and temperature are required for entropy terms')
        if not np.isclose(record['temperature_K'],ts['temperature_K']):
            raise ValueError('Entropy comparison requires equal temperatures')
    result={}
    states={name:{q:r[k] for q,k in [('E','energy_hartree'),('H','H_hartree'),('G','G_hartree')]} for name,r in records.items()}
    states['separated_reactants']={q:states['substrate'][q]+states['Lys'][q] for q in ('E','H','G')}
    for name,v in states.items():
        result[f'{name}_H_hartree']=v['H']
        result[f'{name}_TS_kcal_mol']=(v['H']-v['G'])*KCAL
        result[f'{name}_S_cal_mol_K']=result[f'{name}_TS_kcal_mol']*1000/ts['temperature_K']
    ref=states['separated_reactants']
    for name,label in [('TS','activation'),('product','reaction')]:
        d={q:(states[name][q]-ref[q])*KCAL for q in ('E','H','G')}
        result[f'T_delta_S_{label}_kcal_mol']=d['H']-d['G']
        result[f'minus_T_delta_S_{label}_kcal_mol']=d['G']-d['H']
        result[f'delta_S_{label}_cal_mol_K']=(d['H']-d['G'])*1000/ts['temperature_K']
        result[f'G_minus_E_{label}_kcal_mol']=d['G']-d['E']
        result[f'H_minus_E_{label}_kcal_mol']=d['H']-d['E']
    return result
