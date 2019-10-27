#!/usr/bin/env python3

import subprocess as sp
import os
import os.path
import json
import sys

# Brewer colors
main_colors = ['1f78b4', '33a02c', 'e31a1c', 'ff7f00', '6a3d9a', '8c510a', '01665e', '8e0152', '49006a', '252525', 'aaaa55']
comp_colors = ['a6cee3', 'b2df8a', 'fb9a99', 'fdbf6f', 'cab2d6', 'bf812d', '35978f', 'c51b7d', '7a0177', '525252', 'ffff99']


def main(argv):
    files = []
    for i in range(1,len(argv)):
        files.append(argv[i])

    for fn in files:
        base_name = ".".join(fn.split(".")[0:-1])
        if not os.path.exists(base_name): 
            os.mkdir(base_name)
        
        with open(fn, 'r') as f:
            design_specification = json.load(f)

        domain_lens = {}
        domain_colors = {}
        domain_inds = {}
        domain_colorlist = []
        domain_namelist = []

        j = 0

        for i, domain in enumerate(design_specification['domains']):
            dname = domain['name']
            domain_lens[dname] = len(domain['sequence'])
            if dname[-1] != '*':
                domain_colors[dname] = main_colors[j % len(main_colors)]
                domain_colors[dname + '*'] = comp_colors[j % len(comp_colors)]
                domain_inds[dname] = 2*j
                domain_inds[dname + '*'] = 2*j + 1
                j += 1
                domain_colorlist.append(domain_colors[dname])
                domain_colorlist.append(domain_colors[dname + '*'])
                domain_namelist.append(dname)
                domain_namelist.append(dname + '*')
                # print(dname, main_colors[j-1], comp_colors[j-1])

        strands = {}
        for strand in design_specification['strands']:
            strands[strand['name']] = [cs for cs in strand['domains']]

        structures = {}
        complex_strands = {}
        for structure in design_specification['structures']:
            structures[structure['name']] = structure['structure']
            complex_strands[structure['name']] = [cs for cs in structure['strands']]

        nudraw_cmd = 'drawing-hg/nudraw_python/nudraw.py'

        for name, structure in structures.items():
            full_command = [nudraw_cmd]
            cur_strands = complex_strands[name]
            domain_ind_list = []
            domain_id_list = []
            i = 0
            for strand in cur_strands:
                domains = strands[strand]
                for dom in domains:
                    domain_length = domain_lens[dom]
                    domain_id = domain_inds[dom]
                    domain_ind_list += [str(i)] * domain_length
                    domain_id_list += [str(domain_id)] * domain_length
                    i += 1
            domains_svg = '{}/{}.svg'.format(base_name, name)
            results_svg = '{}/{}_prob.svg'.format(base_name, name)
            
            
            specification_options = [
                    '--drawbaseticks',
                    '--colordomains',
                    '--labeldomains',
                    '--noscale',
                    '--bbwidth=96',
                    '--structure={}'.format(structure),
                    '--domains={}'.format(','.join(domain_id_list)),
                    '--domainids={}'.format(','.join(domain_ind_list)),
                    '--domainnames={}'.format(','.join(domain_namelist)),
                    '--domaincolors={}'.format(','.join(domain_colorlist)),
                    '--svgfile2d={}'.format(domains_svg),
                    ]
            
            results_options = [
                    '--noscale',
                    '--sequence=A',
                    '--colorbaseprob',
                    '--structure={}'.format(structure),
                    '--probfile={}_0_{}.ppairs'.format(base_name, name),
                    '--svgfile2d={}'.format(results_svg),
                    ]
            
            colored_domains = full_command + specification_options
            prob_results = full_command + results_options
            # print(" ".join(colored_domains))
            sp.call(colored_domains)
            # print(" ".join(prob_results))
            sp.call(prob_results)

if __name__ == "__main__":
    main(sys.argv)
