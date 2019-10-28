#!/usr/bin/env python3
# Zhewei Chen
# zchen@caltech.edu
# Edits to nupack drawing script which uses nudraw backend

import subprocess
import os
import os.path
import json
import sys
import matplotlib

# Brewer colors
main_colors = ['1f78b4', '33a02c', 'e31a1c', 'ff7f00', '6a3d9a', '8c510a', '01665e', '8e0152', '49006a', '252525', 'aaaa55']
comp_colors = ['a6cee3', 'b2df8a', 'fb9a99', 'fdbf6f', 'cab2d6', 'bf812d', '35978f', 'c51b7d', '7a0177', '525252', 'ffff99']

def parse_domain_colors(text):
    if text[0]=='[' and text[-1]==']':
        text = text[1:-1].split(',')
    else:
        print('error: syntax for domain colors incorrect')
        help()

    mpl_colors = matplotlib.colors.CSS4_COLORS
    data = {}
    for domain in text:
        x = domain.split(':')
        # map color to ascii value
        if x[1] in mpl_colors.keys():
            data[x[0]] = mpl_colors[x[1]][1:].lower()
        # catch invalid color key
        else:
            print('not valid color name')
            print('the following are valid colors')
            print(mpl_colors.keys())
            help()
    return data

def run_nudraw(fname, colors = {}):
    # get the base name for ppairs file and probability drawing
    base_name = ".".join(fname.split(".")[0:-1])
    if not os.path.exists(base_name): 
        os.mkdir(base_name)
    
    # open the json file and parse it
    with open(fname, 'r') as f:
        design_specification = json.load(f)

    # initial variables
    domain_lens = {}
    domain_colors = {}
    domain_inds = {}
    domain_colorlist = []
    domain_namelist = []

    j = 0

    for i, domain in enumerate(design_specification['domains']):
        dname = domain['name']
        domain_lens[dname] = len(domain['sequence'])
        
        # default if there is custom color for the domain
        if dname[-1] != '*':
            if dname in colors.keys():
                domain_colors[dname] = colors[dname]
                domain_colors[dname + '*'] = colors[dname] # set complement domain as same color
            else:
                domain_colors[dname] = main_colors[j % len(main_colors)]
                domain_colors[dname + '*'] = comp_colors[j % len(comp_colors)]

            domain_inds[dname] = 2*j
            domain_inds[dname + '*'] = 2*j + 1
            j += 1
            domain_colorlist.append(domain_colors[dname])
            domain_colorlist.append(domain_colors[dname + '*'])
            domain_namelist.append(dname)
            domain_namelist.append(dname + '*')

    # parse the strand info
    strands = {}
    for strand in design_specification['strands']:
        strands[strand['name']] = [cs for cs in strand['domains']]

    # parse the structure info
    structures = {}
    complex_strands = {}
    for structure in design_specification['structures']:
        structures[structure['name']] = structure['structure']
        complex_strands[structure['name']] = [cs for cs in structure['strands']]

    # form the commands for nudraw
    nudraw_cmd = '/Users/zchen/Desktop/plab/Computation/nupack_drawing/drawing-hg/nudraw_python/nudraw.py'

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
        subprocess.call(colored_domains)

        prob_results = full_command + results_options
        subprocess.call(prob_results)

def help():
    print('''
-h                          print help
-i <file.json>              input json file from nupack
-c [domain:color]           specify custom color for each domain
                            such as [xRFP:red,xGFP:green]
    ''')
    sys.exit()

def main(argv):
    # show output
    print(argv)

    # initial variables
    colors = ''
    files = []

    for i in range(0,len(argv)):
        if argv[i]=='-h': help()
        elif argv[i] =='-i':
            while(argv[i+1][-5:]=='.json'):
                files.append(argv[i+1])
                i+=1
        elif argv[i] =='-c':
            colors = argv[i+1]
    
    # parse input colors
    colors = parse_domain_colors(colors)
    print('Input files: ',files)
    print('Custom domain colors: ',colors)

    # no input files specified
    if len(files) == 0:
        print('error: no input files given')
        help()
    
    # run if no problems
    for fname in files:
        run_nudraw(fname, colors)

if __name__ == "__main__":
    main(sys.argv)
