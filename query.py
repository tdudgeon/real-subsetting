import sys, time
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs


def gen_fp_morgan2(mol):
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=2048, useChirality=0, useBondTypes=1, useFeatures=0)


def gen_rdk_fp(mol):
    return Chem.RDKFingerprint(mol)


def gen_fp(mol):
    return gen_fp_morgan2(mol)


def read_smiles(file):
    with open(file) as f:
        mol = None
        smi = None
        for line in f:
            if line:
                smi = line.strip()
                mol = Chem.MolFromSmiles(smi)
            yield mol, smi


def read_sdf(file):
    supplr = Chem.SDMolSupplier(file)
    for mol in supplr:
        smi = Chem.MolToSmiles(mol)
        yield mol, smi


def read_file(file):
    if file.endswith('.sdf'):
        return read_sdf(file)
    else:
        return read_smiles(file)


p1 = sys.argv[1]
p2 = sys.argv[2]
threshold1 = float(sys.argv[3])
threshold2 = float(sys.argv[4])
alpha = float(sys.argv[5])

print('Target:', p1, 'Query:', p2, 'MinThreshold:', threshold1, 'MaxThreshold:', threshold2, 'Alpha:', alpha)

t1 = time.time()
fps1 = []
smiles1 = []
for mol, smi in read_file(p1):
    if mol:
        fp = gen_fp(mol)
        fps1.append(fp)
        smiles1.append(smi)

t2 = time.time()
print("Fingerprinting targets took:", (t2 - t1))
print("Number of targets:", len(fps1))

t3 = time.time()
hits = []
best_scores = []
best_target = []
for mol, smi in read_file(p2):
    if mol:
        fp2 = gen_fp(mol)
        scores = DataStructs.BulkTverskySimilarity(fp2, fps1, alpha, 1.0 - alpha)
        best_score = 0.0
        best_target = None

        for score, target in zip(scores, smiles1):
            # print(score)
            if score > best_score:
                best_score = score
                best_target = target
        if threshold1 <= best_score <= threshold2:
            print(smi, best_score)
            print(best_target)
            hits.append(smi)

        # for score in scores:
        #     if threshold1 <= score <= threshold2:
        #         hits.append(smi)
        #         break

t4 = time.time()
print("Queries took:", (t4 - t1))
print('Found', len(hits), 'molecules')

# for hit in hits:
#     print(hit)