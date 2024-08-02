#Reader for TICG simulation data in OVITO.
#
#+-----------------------------------------------+
#| By:                                           |
#|     - Jihun Ahn                               |
#| Date:                                         |
#|     - 2024.07.16                              |
#|                                               |
#+-----------------------------------------------+


from ovito.data import DataCollection, SimulationCell, ParticleType
from ovito.io import FileReaderInterface, import_file
from typing import Callable, Any
import pandas as pd
import numpy as np
import re
import os

def chain_idx(chain_strings):
    with open(chain_strings, 'r') as f:
        for i in range(4):
            f.readline()
        chain_strings = f.readlines()
    chain_idx = []
    idx_num = 0
    for chain in chain_strings:
        N,N_bead,string = chain.split('|')
        for n in range(int(N)):
            chain_idx.extend([idx_num]*int(N_bead))
            idx_num += 1
    return chain_idx


class TICG_Reader(FileReaderInterface):
    def __init__(self):
        self.bond_TICG = None
        self.atom_load = None
    @staticmethod
    def file_extensions():
        return ["ticg"]

    def chain_idx(self,chain_strings):
        with open(chain_strings, 'r') as f:
            for i in range(4):
                f.readline()
            chain_strings = f.readlines()

        chain_idx          = []
        chain_idx_for_bond = []
        idx_num = 0
        for chain in chain_strings:
            N,N_bead,string = chain.split('|')
            for n in range(int(N)):
                chain_idx.extend([idx_num]*int(N_bead))
                chain_idx_for_bond.extend([idx_num]*(int(N_bead)-1))
                idx_num += 1
        return chain_idx, chain_idx_for_bond

    @staticmethod
    def detect(filename: str):
        print(filename)
        try:
            directory = os.path.dirname(filename)
            psf_files = glob.glob(os.path.join(directory, '*.psf'))
            if not psf_files:
                return False
            log_file = os.path.join(directory, 'Simulation.log')
            if not os.path.exists(log_file):
                return False
            ticg_files = glob.glob(os.path.join(directory, '*.ticg'))
            if not ticg_files:
                return False
            return os.path.splitext(filename)[1] == ".ticg"
        except OSError:
            return False


    def scan(self, filename: str, register_frame: Callable[..., None]):
        if os.path.splitext(filename)[1] != ".ticg":
            raise OSError("File is not a TICG file.")
        name = os.path.splitext(filename)[0]
        path = os.path.dirname(filename)
        if not os.path.exists(os.path.join(path, "Simulation.log")):
            print(os.path.join(path, "Simulation.log"))
            raise OSError("Simulation.log file is missing.")
        if not os.path.exists(os.path.join(path, f"{name}.psf")):
            raise OSError("Polymer.psf file is missing.")
        if not os.path.exists(os.path.join(path, "EUV_Input_string")):
            raise OSError("EUV_Input_string file is missing.")
        expr      = re.compile(r'^\d+$')
        label_num = 1
        with open(filename, "r") as f:
            while line := f.readline():
                if re.match(expr, line.strip()): 
                    num_particles = int(line.strip())
                elif "MC simulation of coarse grain block copolymer" in line:
                    offset = f.tell()
                    label = f"Frame {label_num}"
                    label_num += 1
                    register_frame(frame_info=(offset, num_particles), label=label)
        #  print(label_num)
    def parse(self, data: DataCollection, filename: str, frame_info: tuple[int, int], **kwargs: Any):
        starting_byte_offset, num_particles = frame_info
        path = os.path.dirname(filename)
        with open(os.path.join(path, "Simulation.log"), "r") as f:
            log_data = f.readlines()
            for i in range(len(log_data)):
                if 'x:' in log_data[i][:2]:
                    box_x = float(log_data[i].split()[1]) * 10  
                if 'y:' in log_data[i][:2]:
                    box_y = float(log_data[i].split()[1]) * 10
                if 'z:' in log_data[i][:2]:
                    box_z = float(log_data[i].split()[1]) * 10
        cell = data.create_cell([[box_x, 0, 0], [0, box_y, 0], [0, 0, box_z]], pbc=[True, True, False])

        chain_strings = os.path.join(path, "EUV_Input_string")
        chain_idx, chain_idx_for_bond = self.chain_idx(chain_strings)
        #  print(f.seek(starting_byte_offset))
        name_dict = {
                0:  ['protection'   , 0.25, [254., 85., 85.]] ,
                1:  ['deprotection' , 0.25, [85., 85., 254.]] ,
                2:  ['backbone'     , 0.25, [254., 248., 172.]] ,
                3:  ['developer'    , 0.1 , [201., 200., 200.]] ,
                4:  ['PAG'          , 0.25, [156., 254., 139.]] ,
                5:  ['Quen'         , 0.25, [254., 93., 241.]] ,
                6:  ['Acid'         , 0.25, [0., 255., 255.]] ,
                     }
        with open(filename, "r") as f:
            f.seek(starting_byte_offset)
            df = pd.read_csv(f, sep='\s+', header=None, nrows=num_particles,names=['type', 'x', 'y', 'z'])
            df['chain_idx'] = 0
            df['dev'] = 0
            df.iloc[:len(chain_idx),4] = chain_idx
            #  df['type'] = df['type'].replace(name_dict)
            df.iloc[len(chain_idx):,4] = -df.iloc[len(chain_idx):,0]
            df.loc[df['type'] == 3, 'dev'] = 1
            particles     = data.create_particles(count=num_particles,vis_params={'radius': 0.5})
            particle_type = particles.create_property("Particle Type", data=df['type'].values)
            positions     = particles.create_property("Position", data=df[['x', 'y', 'z']].values)
            chain_type    = particles.create_property("Chain Type", data=df['chain_idx'].values)
            developer     = particles.create_property("Developer", data=df['dev'].values)
        if self.atom_load is None:
            for particle_id, args in name_dict.items():
                #  particle_type.add_particle_type(name)
                name   = args[0]
                radius = args[1]
                #  color  = args[2]
                color = np.array(args[2])/255

                particle_type.types.append(ParticleType(id = particle_id, name = name, radius = radius, color = color))
                self.atom_load = True
        #check if the bond information is available

        if self.bond_TICG is not None:
            bonds = self.bond_TICG
            bond = particles.create_bonds(vis_params={'width': 0.1})
            bond.create_property("Topology", data=bonds)
            bond.create_property("Type", data=chain_idx_for_bond)
            print("bond information is available")
            print("Load complete")
        else:
            print("bond information is not available")
            print("Calculating bond information")
            with open(os.path.join(path, f"{os.path.splitext(filename)[0]}.psf"), "r") as f:
                psf_data = f.readlines()
                #  print(psf_data)
                start_line = num_particles + 9
                for i in range(len(chain_idx)):
                    if '!NTHETA' in psf_data[-i]:
                        end_line = len(psf_data) - i
                        break
                bonds = psf_data[start_line:end_line]
                bonds = pd.DataFrame([i.split() for i in bonds], columns=['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2'])
                bonds_list = bonds[['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2']].values.reshape(-1, 2)
                bonds = pd.DataFrame(bonds_list, columns=['idx1', 'idx2']).dropna().astype(int).values - 1
                self.bond_TICG = bonds
                bond = particles.create_bonds(vis_params={'width': 0.1})
                bond.create_property("Topology", data=bonds)
                bond.create_property("Type", data=chain_idx_for_bond)
                print("Load complete")
        

if __name__ == "__main__":
    pipeline = import_file("../../EUV/Development/N64/bead_out_develop.ticg", input_format=TICG_Reader)
    string = chain_idx("../../EUV/Development/N64/EUV_Input_string")
    #  print(len(string))
    #  for frame in range(pipeline.source.num_frames):
    #      data = pipeline.compute(frame)

