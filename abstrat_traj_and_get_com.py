# -*- coding:utf-8 -*-
# Author ：Xinyue Sun
# Time : 2022/11/01
# Description: abstract frames what you want from trajectory files and obtain center of mass of ions (cations and anions) from abstracted trajectory files(.xyz)

import numpy as np
import pandas as pd
import MDAnalysis as mda
import os, re

## abstract frames what you want
class read_trj():
    def __init__(self, trj_file):
        self.trj_file = trj_file

    def read_trj(self):
        frames = 0
        timesteps = []

        if os.path.exists(self.trj_file):
            f2 = open(self.trj_file, 'r')
            sum_of_atoms = f2.readline()
            sum_of_atoms = int(sum_of_atoms)
            for line in f2:
                if "Atoms." in line:
                    frames += 1
                    timesteps.append(line)
            f2.close()
            return frames, timesteps
        else:
            raise IOError("the trjactory file is not existed!")

    def abstract_trj(self, frames):
        sum_of_atoms = 0
        timesteps = []
        f3 = open(self.trj_file, 'r')
        sum_of_atoms = f3.readline()
        sum_of_atoms = int(sum_of_atoms)
        f3.seek(0)
        f4 = open("abstracted_traj", 'a+')
        frames = [int(x) for x in frames]
        print(frames)
        f, timesteps = self.read_trj()

        for line in f3:
            for i in frames:
                timestep = timesteps[i - 1]
                # print(timestep)
                if timestep in line:
                    f4.write(str(sum_of_atoms) + '\n')
                    f4.write(timestep)
                    for i in range(sum_of_atoms):
                        line = f3.readline()
                        f4.write(line)


# read molecular information from mmol files
if os.path.exists('abstracted_traj'):
    os.remove("./abstracted_traj")
trj_file = input('please input the trj_file name:')
read_trj_inst = read_trj(trj_file)
frames, steps = read_trj_inst.read_trj()
frames = int(frames)
print(frames)

abstracted_frames = input(
    "please input the frames you want to choose: all, one or range such as all or 1 2 3, or 1-5: ")
if abstracted_frames == "all":
    new_list = np.arange(1, frames + 1, 1)
    new_list.tolist()
    print(new_list)
    read_trj_inst.abstract_trj(new_list)

elif "-" in abstracted_frames:
    list = abstracted_frames.split("-")
    new_list = np.arange(int(list[0]), int(list[1]) + 1, 1)
    new_list = new_list.tolist()
    print(new_list)
    read_trj_inst.abstract_trj(new_list)

else:
    new_list = abstracted_frames.split(' ')
    print(new_list)
    read_trj_inst.abstract_trj(new_list)


## get center of mass
# Relative atomic mass of the element
atomic_wt = {'H': 1.008, 'Li': 6.941, 'B': 10.811, 'C': 12.011,
             'N': 14.007, 'O': 15.9994, 'F': 18.998, 'Ne': 20.180,
             'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086,
             'P': 30.974, 'S': 32, 'Cl': 35.453, 'Ar': 39.948,
             'K': 39.098, 'Ca': 40.078, 'Ti': 47.867, 'Fe': 55.845,
             'Zn': 65.38, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
             'Mo': 95.96, 'Ru': 101.07, 'Sn': 118.710, 'Te': 127.60,
             'I': 126.904, 'Xe': 131.293}

def fetchList(filePath, splitString, indexingStart, indexingEnd='None', skipBlankSplit=False, *args):
    readMode = False
    openedFile = open(filePath, 'r')
    listArray = []
    currentLine = ''
    for line in openedFile:
        previousLine = currentLine
        currentLine = str(line)
        if currentLine.isspace() and readMode is True:  # In each line, check if it's purely white space
            # If the previous line doesn't countain the spliting string (i.e. this isn't the line directly after the split), stop reading; only occurs if the skipBlankSplit option is set to True
            if splitString not in previousLine or skipBlankSplit is False:
                readMode = False
        if currentLine[0] == '[':  # The same is true if it's the start of a new section, signalled by '['
            readMode = False
        if readMode is True:
            if line[0] != ';' and currentLine.isspace() is False:  # If the line isn't a comment
                if indexingEnd == 'None':
                    # Just split with a start point, and take the rest
                    importantLine = (currentLine.split())[indexingStart:]
                else:
                    # Define the bigt we care about as the area between index start and end
                    importantLine = (currentLine.split())[indexingStart:indexingEnd]
                # If any additional arguments were given (and they're numbers), also add this specific index in
                for ar in args:
                    if is_number(ar):
                        importantLine += [(currentLine.split())[ar]]
                listArray.append(importantLine)  # Append the important line to the list array
        if splitString in currentLine:  # Check at the end because we only want to start splitting after the string; if found, it'll split all further lines till it gets to white space
            readMode = True
    openedFile.close()
    return listArray  # Return the array containing all associations

def read_traj(topo, traj, ilnum, waternum):
    u = mda.Universe(topo, atom_style='id resid type charge x y z')
    # add elements
    elements = []
    for i in range(len(u.atoms)):
        elements.append(chemical_symbols[np.where(types == int(u.atoms[i].type))[0][0]])
    u.add_TopologyAttr('element', values=elements)

    # add mol. name
    if waternum == 0:
        resnames = ['cation'] * ilnum + ['anion'] * ilnum
    else:
        resnames = ['cation'] * ilnum + ['anion'] * ilnum + ['water'] * waternum
    u.add_TopologyAttr('resnames', values=resnames)

    u.load_new(traj, format="LAMMPSDUMP", timeunit="fs", dt=10000)
    workflow = [transformations.unwrap(u.atoms)]
    u.trajectory.add_transformations(*workflow)
    return u


def atomic_symbol(name):
    if name[:2] in atomic_wt:
        return name[:2]
    elif name[0] in atomic_wt:
        return name[0]
    else:
        print('warning: unknown symbol for atom '  +name)
        return name


directory = './'
filename = 'data.lmp'
## elements
elements = fetchList(filename, 'Masses', 0,4, skipBlankSplit=True)
digitspattern = r'#'
types = []
chemical_symbols = []
for element in elements:
    types.append(int(element[0]))
    txt = re.sub(digitspattern, '', element[-1])
    chemical_symbols.append(atomic_symbol(txt))
types = np.array(types)
print('types:            ',types)
print('chemical_symbols: ',chemical_symbols)

il_num=500
# full: atom-ID molecule-ID atom-type q x y z
u = mda.Universe('result.data', atom_style='id resid type charge x y z')
# resnames
resnames = ['cati']*il_num+['anio']*il_num+['MXene']*1+['MXene']*1
u.add_TopologyAttr('resnames',values=resnames)
u.load_new("abstracted_traj", format="XYZ",timeunit="fs",dt=10000)
print('dimensions:',u.dimensions)
print('frames',u.trajectory.n_frames)

cations = u.select_atoms('resname cati') # select all cations
resids_ca = np.unique(cations.atoms.resids) # extract resid of all cations
cation_coms = []
for ts in u.trajectory:
    for i in resids_ca:
        sel = u.select_atoms('resid %d'%i)
        com = sel.center_of_mass()
        cation_coms.append([com[0],com[1],com[2]])
cation_com_pd = pd.DataFrame(cation_coms)
with open('cacom','w') as f:
    f.write('\n')
    np.savetxt(f,cation_com_pd.values)


anions = u.select_atoms('resname anio') # select all anions
resids_an = np.unique(anions.atoms.resids) # extract resid of all anions
anion_coms = []
for ts in u.trajectory:
    for i in resids_an:
        sel = u.select_atoms('resid %d'%i)
        com = sel.center_of_mass()
        anion_coms.append([com[0],com[1],com[2]])
anion_com_pd = pd.DataFrame(anion_coms)
with open('ancom','w') as f:
    f.write('\n')
    np.savetxt(f,anion_com_pd.values)





