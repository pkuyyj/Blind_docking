import os
import json
import csv
import glob
from rmsdfn import rmsd, fast_rmsd, get_coor_from_sdf

from argparse import FileType, ArgumentParser
parser = ArgumentParser()
parser.add_argument('--input_path', type=str, default='input', help='Path to folder with testing ligand-protein structures')
parser.add_argument('--diffdock_result_path', type=str, default='diffdock/results0129', help='Path to folder with testing ligand-protein structures')

parser.add_argument('--results_path', type=str, default='results', help='Path to save results')

parser.add_argument('--seed', type=str, default='666', help='Random seed for docking')
parser.add_argument('--date', type=str, default='0217', help='date for result folder naming')
parser.add_argument('--search_mode', type=str, default='32', help='exhaustiveness for docking searching')
parser.add_argument('--ori_pocket', action='store_true', default=False, help='use exact pocket of Diffdock output without normalize')
parser.add_argument('--pocket_size', type=float, default=30, help='normalized pocket size in Angstrom')

args = parser.parse_args()

input_path = args.input_path

seed = args.seed
date = args.date
sf = 'vina'
# mode = "detail"
mode = args.search_mode
normal_pocket = not args.ori_pocket
pocket_size = args.pocket_size
if normal_pocket:
    output_path = os.path.join(args.results_path, mode+'_'+sf+"_diffdock_box_pdbqt_vina_seed"+str(seed) +"_normal"+str(pocket_size))+'_date'+date
else:
    output_path = os.path.join(args.results_path, mode+'_'+sf+"_diffdock_box_pdbqt_vina_seed"+str(seed))+'_date'+date

os.system("mkdir "+output_path)

out_csv = os.path.join(args.results_path, "rmsd_results_{}_{}_diffdock_box_pdbqt_vina_seed{}_date{}.csv".format(mode,sf,seed,date))
with open(out_csv,"w") as f:
    f.write("pdb,rmsd\n")
    
configs = sorted(glob.glob(os.path.join(args.diffdock_result_path,"*/rank1.sdf")),reverse=True)


for config in configs:
    pdb = config.split('/')[-2].split('-')[-2]
    with open(config, "r") as f:
        lines = f.readlines()
    
    ligand = os.path.join(input_path,'{}/{}_ligand.pdbqt'.format(pdb,pdb))    
    
    ori_coords, _ = get_coor_from_sdf(config)
    
    center_max = [-1e4,-1e4,-1e4]
    center_min = [1e4,1e4,1e4]
    center = [0,0,0]
    size = [0,0,0]
    for co in ori_coords:
        for i in range(3):
            center_max[i] = max(center_max[i], co[i])
            center_min[i] = min(center_min[i], co[i])
    
    for i in range(3):
        center[i] = (center_max[i] + center_min[i]) / 2
        size[i] = (center_max[i] - center_min[i])  + 5
    
    
    new_config = os.path.join(output_path,pdb+"_config.txt")
    with open(new_config,"w") as f:
        f.write("center_x = {}\ncenter_y = {}\ncenter_z = {}\n".format(center[0],center[1],center[2]))
        if normal_pocket:
            f.write("size_x = {}\nsize_y = {}\nsize_z = {}\n".format(pocket_size,pocket_size,pocket_size))
        else:
            f.write("size_x = {}\nsize_y = {}\nsize_z = {}\n".format(min(size[0],60),min(size[1],60),min(size[2],60)))
        if mode in ['detail', 'balance']:
            f.write("search_mode = {}\n".format(mode))
        else:
            f.write("exhaustiveness = {}\n".format(mode.split('_')[0]))
        
        f.write("receptor = {}/{}_protein.pdbqt\nligand = {}\n".format(os.path.join(input_path,pdb),pdb,ligand))
        f.write("out = {}\nseed = {}\nnum_modes = 9\nverbosity = 1\n".format(os.path.join(output_path,pdb+"_ligand_out.pdbqt"), seed))
        
    if not os.path.isfile(os.path.join(output_path,pdb+"_ligand_out.pdbqt")):
        os.system("vina --config "+new_config)
    if not os.path.isfile(os.path.join(output_path,pdb+"_ligand_out.pdbqt")):
        print("error in "+ pdb)
        with open(out_csv, "a") as out_f:
            out_f.write(pdb+",999\n")
        continue
    rmsd_ = fast_rmsd(ligand, os.path.join(output_path, pdb+"_ligand_out.pdbqt"))
    with open(out_csv, "a") as out_f:
        out_f.write(pdb+','+str(rmsd_)+"\n")
        
