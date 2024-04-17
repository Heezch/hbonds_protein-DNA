import numpy as np
import mdtraj as md

def check_hbond(top, hbond, protein_chainid=0):
    donor_id = top.atom(hbond[0]).residue.chain.index
    acceptor_id = top.atom(hbond[2]).residue.chain.index

    condition_protein_donor = (donor_id == protein_chainid) and (acceptor_id != protein_chainid)
    condition_protein_acceptor = (donor_id != protein_chainid) and (acceptor_id == protein_chainid)

    return condition_protein_donor or condition_protein_acceptor

def filter_hbonds_for_frame(hbonds, top, protein_chainid=0):
    return [hbond for hbond in hbonds if check_hbond(top, hbond, protein_chainid=protein_chainid)]

def load_trajectory_and_compute_hbonds(file_path):
    traj = md.load(file_path).remove_solvent()
    hbond_traj = md.wernet_nilsson(traj)  # returns a list of tuples of hydrogen bond atom indices (donor, hydrogen, acceptor)
    top = traj.top
    return traj, hbond_traj, top

def process_hbonds(traj, hbond_traj, top, protein_chainid=2):
    filtered_hbonds = [filter_hbonds_for_frame(hbonds, top, protein_chainid) for hbonds in hbond_traj]
    filtered_hbonds_array = np.array(filtered_hbonds, dtype=object)
    return filtered_hbonds_array

def print_hbonds(traj, filtered_hbonds_array):
    print('List of hydrogen bonds in the first frame:')
    atoms = traj.top._atoms
    for h in filtered_hbonds_array[0]:
        print(atoms[h[0]], atoms[h[1]], atoms[h[2]])

    n_hbonds = [len(hbonds) for hbonds in filtered_hbonds_array]
    print(f'There are {n_hbonds[0]} hydrogen bonds in the first frame between the protein and DNA')

def main():
    # Specify your trajectory file path
    file_path = './your_protein_DNA.pdb'
    traj, hbond_traj, top = load_trajectory_and_compute_hbonds(file_path)
    
    # Print the chains and residues
    for c in top.chains:
        r = [r for r in c.residues]
        print(c.index, r)
    
    filtered_hbonds_array = process_hbonds(traj, hbond_traj, top, protein_chainid=2)
    print_hbonds(traj, filtered_hbonds_array)

if __name__ == '__main__':
    main()
