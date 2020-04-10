from monty.serialization import loadfn, dumpfn
from pymatgen import Structure, Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.util.plotting import periodic_table_heatmap
from itertools import groupby

scan_entries = loadfn('element_scan_entries_2020-03-16.json')
gga_entries = loadfn('element_gga_entries_2020-03-10.json')

def sort_by_element(entries_list):
    elem_dict = dict()
    for entry in entries_list:
        element = entry.structure.composition.reduced_formula
        if not element in elem_dict.keys():
            elem_dict[element] = list([entry])
        else:
            elem_dict[element].append(entry)
    return elem_dict

scan_sorted, gga_sorted = sort_by_element(scan_entries), sort_by_element(gga_entries)

missing_from_scan = list()
for elem in gga_sorted.keys():
    if elem not in scan_sorted.keys():
        missing_from_scan.append(elem)

def get_ground_state(entries_dict):
    gs_dict = dict()
    for element in entries_dict.keys():
        entry_set = EntrySet(entries_dict[element])
        entry_set.remove_non_ground_states()
        gs_dict[element] = list(entry_set)[0]
    return gs_dict

scan_gs_dict, gga_gs_dict = get_ground_state(scan_sorted), get_ground_state(gga_sorted)

def get_diff_dict(scan_gs_dict, gga_gs_dict, ltol=0.2, stol=0.3, angle_tol=5, allow_subset=False):
    matcher = StructureMatcher(ltol, stol, angle_tol, allow_subset)
    diff_dict = dict()
    for element in scan_gs_dict.keys():
        scan_struct, gga_struct = scan_gs_dict[element].structure, gga_gs_dict[element].structure
        if not matcher.fit(scan_struct, gga_struct):
            diff_dict[element] = (scan_gs_dict[element], gga_gs_dict[element])
    return diff_dict

l, s, a, = 0.2, 0.3, 5
tolerances = [{'ltol':l  , 'stol':s  , 'angle_tol':a  , 'allow_subset':False},
              {'ltol':2*l, 'stol':s  , 'angle_tol':a  , 'allow_subset':False},
              {'ltol':l  , 'stol':2*s, 'angle_tol':a  , 'allow_subset':False},
              {'ltol':l  , 'stol':s  , 'angle_tol':a*2, 'allow_subset':False},
              {'ltol':2*l, 'stol':s  , 'angle_tol':a*2, 'allow_subset':False},
              {'ltol':2*l, 'stol':2*s, 'angle_tol':a  , 'allow_subset':False},
              {'ltol':l  , 'stol':2*s, 'angle_tol':a*2, 'allow_subset':False},
              {'ltol':2*l, 'stol':2*s, 'angle_tol':a*2, 'allow_subset':False},
              {'ltol':2*l, 'stol':2*s, 'angle_tol':a*2, 'allow_subset':True}]

entry_configs = list()
for config in tolerances:
    diff = get_diff_dict(scan_gs_dict, gga_gs_dict, config['ltol'], config['stol'], config['angle_tol'], config['allow_subset'])
    entry_configs.append(diff)

energy_configs = dict()
# for diff in entry_configs:
#     d = dict()
#     for elem in diff.keys():
#         d[elem] = diff[elem][0].energy - diff[elem][1].energy
#     energy_configs.append(d)
diff = entry_configs[-1]
for elem in diff.keys():
    energy_configs[elem] = diff[elem][0].energy - diff[elem][1].energy

periodic_table_heatmap(energy_configs, show_plot=True, cmap='Spectral')


# blank_color='white', cbar_label=str(config[-1]),
# diffs = list()

# diff0 = get_diff_list(scan_gs_dict, gga_gs_dict)
# diff1 = get_diff_list(scan_gs_dict, gga_gs_dict, ltol=0.4)
# diff2 = get_diff_list(scan_gs_dict, gga_gs_dict, stol=0.6)
# diff3 = get_diff_list(scan_gs_dict, gga_gs_dict, angle_tol=10)
# diff4 = get_diff_list(scan_gs_dict, gga_gs_dict, ltol=0.4, angle_tol=10)
# diff5 = get_diff_list(scan_gs_dict, gga_gs_dict, ltol=0.4, stol=0.6)
# diff6 = get_diff_list(scan_gs_dict, gga_gs_dict, stol=0.6, angle_tol=10)
# diff7 = get_diff_list(scan_gs_dict, gga_gs_dict, ltol=0.4, stol=0.6, angle_tol=10)
# diff8 = get_diff_list(scan_gs_dict, gga_gs_dict, ltol=0.4, stol=0.6, angle_tol=10, allow_subset=True)

# l1 = [' ' if elem in diff1[0] else 'X' for elem in diff0[0]]
# l2 = [' ' if elem in diff2[0] else 'X' for elem in diff0[0]]
# l3 = [' ' if elem in diff3[0] else 'X' for elem in diff0[0]]
# l4 = [' ' if elem in diff4[0] else 'X' for elem in diff0[0]]
# l5 = [' ' if elem in diff5[0] else 'X' for elem in diff0[0]]
# l6 = [' ' if elem in diff6[0] else 'X' for elem in diff0[0]]
# l7 = [' ' if elem in diff7[0] else 'X' for elem in diff0[0]]
# l8 = [' ' if elem in diff8[0] else 'X' for elem in diff0[0]]


def write_to_cif(scan_structs, gga_structs):
    for struct in scan_structs:
        struct.to(fmt=".cif", filename='cif_files/scan/{0}_scan.cif'.format(struct.composition.reduced_formula))
    for struct in gga_structs:
        struct.to(fmt=".cif", filename='cif_files/gga/{0}_gga.cif'.format(struct.composition.reduced_formula))

write_to_cif([c.structure for c in scan_gs_dict.values()], [c.structure for c in gga_gs_dict.values()])