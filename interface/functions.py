import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import geopandas
import descartes
import ipywidgets as widgets
import matplotlib.gridspec as gridspec
from collections import Counter
import re
import webbrowser
import time
import pickle
import math
from IPython.display import display
import seaborn as sns
sns.set_style('white')

pd.set_option('display.max_rows', 1000)

def sppnames():
    sppnames = {'cow':'bos_taurus',
     'dog': 'canis_lupus',
     'goat': 'capra_hircus',
     'horse': 'equus_caballus',
     'cat':'felis_catus',
     'rabbit': 'oryctolagus_cuniculus',
     'sheep': 'ovis_aries',
     'pig': 'sus_scrofa'}
    return (sppnames)


def FastaToDict(infile, lengths_only=False, spliton=False):
    D = dict()
    seq = []
    nam = ""
    x = 0
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) != 0:
                if line[0] == ">":
                    if len(seq) != 0:
                        seq = "".join(seq)
                        if lengths_only:
                            D[nam] = len(seq)
                        else:
                            D[nam] = seq
                        seq = []
                    nam = line.replace(">", "")
                    if spliton is not False:
                        nam = nam.split(spliton)[0]
                else:
                    seq.append(line)
            x += 1
    seq = "".join(seq)
    if x != 0:
        D[nam] = seq
    return D


def getRefLengths():
    lengths = {'Feline_leukemia_virus': 8448,
               'Mycoplasma_haemofelis': 1155937,
               'Bartonella_henselae': 1931047,
               'Borna_disease_virus': 8910,
               'Equine_infectious_anemia_virus': 8359,
               'Neorickettsia_risticii': 879977,
               'Bartonella_vinsonii': 1995085,
               'Canine_parvovirus': 5323,
               'Rabies_lyssavirus': 11932,
               'Brucella_abortus': 3244234,
               'Bovine_leukemia_virus': 8419,
               'Foot_and_mouth_disease_virus': 8201,
               'Anaplasma_phagocytophilum': 1471282,
               'Escherichia_coli': 4648643,
               'Visna_virus': 9202,
               'Classical_swine_fever_virus': 12104,
               'Hepatitis_E_virus': 7223,
               'Brucella_suis': 3306334,
               'Influenza_A_virus': 13588,
               'Myxoma_virus': 161773,
               'Porcine_reproductive_and_respiratory_syndrome_virus': 15411,
               'Torque_teno_sus_virus': 2878}
    return (lengths)

def getData():
    data = pd.read_csv("../data/clean_data_geo_host.tsv", sep="\t")
    # data = pd.read_csv("https://raw.githubusercontent.com/KatyBrown/ibl_software/master/data/clean_data_geo_host.tsv", sep="\t")
    return (data)

def getSubtab(microbes=[], hosts=[], typs=[], countries=[]):
    if not isinstance(microbes, list):
        microbes = [microbes]
    if not isinstance(hosts, list):
        hosts = [hosts]
    if not isinstance(typs, list):
        typs = [typs]
    if not isinstance(countries, list):
        countries = [countries]
    data = getData()
    if len(microbes) == 0:
        microbes = list(set(data['microbe_name']))
    if len(hosts) == 0:
        hosts = list(set(data['host_name']))
    if len(typs) == 0:
        typs = list(set(data['microbe_type']))
    if len(countries) == 0:
        countries = list(set(data['country']))
    subtab = data[data['microbe_name'].isin(microbes) & (data['host_name'].isin(hosts)) & data['microbe_type'].isin(typs) & data['country'].isin(countries)]
    return (subtab)

def getPathNames():
    F = FastaToDict("../data/references.fasta")
    m = list(F.keys())
    m = [a.lower() for a in m]
    return (m)

def getOrder():
    return (['cat', 'horse', 'dog', 'cow', 'sheep', 'pig'])

def getSppAnswers():
    spp_answers = {'1': 'cat',
                   '2': 'horse',
                   '3': 'dog',
                   '4': 'cow',
                   '5': 'sheep',
                   '6': 'pig'}
    return (spp_answers)

def reverseAnswers():
    spp_answers = {'cat': '1',
                   'horse': '2',
                   'dog': '3',
                   'cow': '4',
                   'sheep': '5',
                   'pig': '6'}
    return (spp_answers) 
  
def getPathogens():
    pathogens = {'1':['feline_leukemia_virus',
                      'mycoplasma_haemofelis',
                      'bartonella_henselae'],
                 '2': ['borna_disease_virus',
                       'equine_infectious_anemia_virus',
                       'neorickettsia_risticii'],
                 '3': ['bartonella_vinsonii',
                       'canine_parvovirus',
                       'rabies_lyssavirus'],
                 '4': ['brucella_abortus',
                       'bovine_leukemia_virus',
                       'foot_and_mouth_disease_virus'],
                 '5': ['anaplasma_phagocytophilum',
                       'escherichia_coli',
                       'visna_virus'],
                 '6': ['classical_swine_fever_virus',
                       'hepatitis_e_virus',
                       'brucella_suis']}
    return (pathogens)
def getSppLoc():
    spp_loc = {'cat': 'brazil',
               'horse': 'usa',
                'dog': 'thailand',
                'cow': 'turkey',
                'sheep': 'china',
                'pig': 'india'}
    return (spp_loc)

def getLoc():
    loc = {'1': 'brazil',
           '2': 'united_states_of_america',
           '3': 'thailand',
           '4': 'turkey',
           '5': 'china',
           '6': 'india'}
    return (loc)

def getDatasetReads(dataset, typ='path'):
    if typ == 'path':
        F = FastaToDict("../data/sample_%s_reads.fasta" % dataset)
    elif typ == 'food':
        F = FastaToDict("../data/sample_%s_reads_food.fasta" % dataset)
    elif typ == 'digest':
        F = FastaToDict("../data/sample_%s_reads_digest.fasta" % dataset)
    elif typ == 'env':
        F = FastaToDict("../data/sample_%s_reads_env.fasta" % dataset)
    return (F)
    
def getColumn(column, hosts=[], typs=[], microbes=[], countries=[]):
    s = getSubtab(hosts=hosts, typs=typs, microbes=microbes, countries=countries)
    return (set(s[column]))

def getDropdown(column, hosts=[], typs=[], microbes=[], countries=[], plus=[]):
    if not isinstance(plus, list):
        plus = [plus]
    opts = ['select option'] + plus + sorted(list(getColumn(column, hosts=hosts, typs=typs, microbes=microbes, countries=countries)))
    opts = [o.replace("_", " ").capitalize().strip() for o in opts]
    w = widgets.Dropdown(options=opts)
    return (w)
    
def plotPathDist(microbe, host, save=False):
    microben = microbe
    hostn = host
    microbe = microbe.replace(" ", "_").lower()
    host = host.replace(" ", "_").lower()
    if host == 'select_option' or microbe == 'select_option':
        return None
    if host == 'all_hosts':
        host = []
    if microbe == 'all_microbes':
        microbe = []
   # world = pd.read_csv("world.tsv", sep="\t")
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world['name_lower'] = world['name'].apply(lambda x: x.lower())
    world = world[['name_lower', 'geometry']]
    f = plt.figure(figsize=(20, 10), dpi=200)
    a = f.add_subplot('111')
    a.set_facecolor('#B5F0FC')
    a.set_xlim(-180, 180)
    a.set_ylim(-90, 90)
    a.text(-180, 95, '%s, %s' % (microben, hostn), fontsize=20)
    subtab = getSubtab(hosts=host, microbes=microbe)
    if microbe == [] and host == []:
        s2 = subtab.groupby(['country']).sum()
        s3 = s2[['number_of_samples']]
        s3.index = range(len(s3))
        s3['country'] = s2.index
        subtab = s3        
    elif microbe == []:
        s2 = subtab.groupby(['host_name', 'country']).sum()
        s3 = s2[['number_of_samples']]
        s3.index = range(len(s3))
        s3['host_name'] = [a for a, b in s2.index]
        s3['country'] = [b for a, b in s2.index]
        subtab = s3
    elif host == []:
        s2 = subtab.groupby(['microbe_name', 'country']).sum()
        s3 = s2[['number_of_samples']]
        s3.index = range(len(s3))
        s3['microbe_name'] = [a for a, b in s2.index]
        s3['country'] = [b for a, b in s2.index]
        subtab = s3
    if len(subtab) != 0:
        colours = matplotlib.cm.Reds(np.linspace(0, 1, max(subtab['number_of_samples'])+1))
        colours = [matplotlib.colors.rgb2hex(x) for x in colours]
    else:
        colours = []
    s = []
    for ind in world.index.values:
        row = world.loc[ind]
        country = row['name_lower']
        p = row['geometry']
        if len(subtab) == 0:
            a.add_patch(descartes.PolygonPatch(p, fc='white', ec='black'))
            a.text(0, 0, "No Results!", color='red', fontsize=40,
                   ha='center', va='center')
        else:
            subsubtab = subtab[subtab['country'] == country]
            if len(subsubtab) != 0:
                row2 = subsubtab.loc[subsubtab.index.values[0]]
                count = row2['number_of_samples']
                colour = colours[int(row2['number_of_samples'])]
                if p.area < 5:
                    c = p.representative_point()
                    a.scatter(c.x, c.y, s=200, facecolor=colour, edgecolor='black', marker="*", zorder=10)
                s.append([country, count])
                a.add_patch(descartes.PolygonPatch(p, fc=colour, ec='black'))
            else:
                a.add_patch(descartes.PolygonPatch(p, fc='white', ec='black'))
    a.set_axis_off()
    if save:
        f.savefig("../interface/%s_%s.png" % (microbe, host), dpi=200, bbox_inches='tight')
    df = pd.DataFrame(s, columns=['Country', 'Number of Samples'])
    df['Country'] = df['Country'].apply(lambda x: x.title())
    df = df.sort_values('Country')
    df['Number of Samples'] = df['Number of Samples'].astype(int)
    return(df)


def splitN(n):
    n = str(n)
    pi = len(n)
    x = []
    for i in range(len(n), 0, -3):
        x.append(n[i:pi])
        pi = i
    x.append(n[:pi])
    return (",".join(x[::-1][:-1]))


def sampleReads(data, n):
    s = np.random.choice(list(data.keys()), n)
    for key in s:
        print (">%s\n%s" % (key, data[key]))

def goToReads(dataset):
    webbrowser.open("../data/sample_%s_reads.fasta" % dataset)

def mapReadsDisplay(dataset, f, typ='path', diet=None):
    fasta = getDatasetReads(dataset, typ)
    refs = FastaToDict(f)
    print ("Loading data....\n")
    time.sleep(3)
    print ("Loaded 40,000 reads")
    time.sleep(3)
    print ("Mapping reads...\n")
    time.sleep(3)
    starts = dict()
    s = list(refs.keys())
    np.random.shuffle(s)
    t = 0
    if diet:
        f = pd.read_csv("../data/food_typs.tsv", sep="\t", header=None)
        if diet == 'carnivore':
            f = f[f[1] == 'c']
        elif diet == 'herbivore':
            f = f[f[1] == 'h']
        f[0] = [x.capitalize().replace(" ", "_") for x in f[0]]
        s = list(f[0])
    for i, nam in enumerate(s):
        nam2 = nam.replace("_", " ")
        print ("%i / %i: Mapping reads to %s" % (i+1, len(s), nam2))
        time.sleep(1)
        ref = refs[nam]
        starts.setdefault(nam, [])
        for seq in fasta.values():
            pos = ref.find(seq)
            if pos != -1:
                starts[nam].append(pos)
        print ("%i reads identified matching %s" % (len(starts[nam]), nam2))
        time.sleep(1)
        t += len(starts[nam])
        print ("%i reads identified in total\n" % (t))
        time.sleep(1)
    print ("Identified %i reads!" % t)
    time.sleep(1)
    print ("Done")

def mapReadsActual(dataset, typ='path', diet=None):
    if typ == 'path':
        f = 'references.fasta'
        p = '.pickle'
    elif typ == 'digest':
        f = 'digest.fasta'
        p = '_digest.pickle'
    elif typ == 'env':
        f = 'env.fasta'
        p = '_env.fasta'
    elif typ == 'food':
        f = 'food.fasta'
        p = '_food.fasta'  
    refs = FastaToDict("../data/%s" % f)
    fasta = getDatasetReads(dataset, typ=typ)
    starts = dict()
    s = list(refs.keys())
    if diet:
        f = pd.read_csv("../data/food_typs.tsv", sep="\t", header=None)
        if diet == 'carnivore':
            f = f[f[1] == 'c']
        elif diet == 'herbivore':
            f = f[f[1] == 'h']
        s = list(f[0])
    np.random.shuffle(s)
    for i, nam in enumerate(s):
        ref = refs[nam]
        starts.setdefault(nam, [])
        for seq in fasta.values():
            pos = ref.find(seq)
            if pos != -1:
                starts[nam].append(pos)
    return (starts)
    
def showMapping(dataset, typ='path', diet=None, save=False):
    if typ == 'path':
        f = 'references.fasta'
        p = '.pickle'
    elif typ == 'digest':
        f = 'digest.fasta'
        p = '_digest.pickle'
    elif typ == 'env':
        f = 'env.fasta'
        p = '_env.fasta'
    elif typ == 'food':
        f = 'food.fasta'
        p = '_food.fasta'
    refs = FastaToDict("../data/%s" % f)
    o = open("../data/%s%s" % (dataset, p), "rb")
    r = pickle.load(o)
    o.close()
    results = dict()
    refnames = list(refs.keys())
    if diet:
        f = pd.read_csv("../data/food_typs.tsv", sep="\t", header=None)
        if diet == 'carnivore':
            f = f[f[1] == 'c']
        elif diet == 'herbivore':
            f = f[f[1] == 'h']
        refnames = list(f[0])
    np.random.shuffle(refnames)
    for ref in refnames:
        if len(r[ref]) != 0:
            L = r[ref]
            x = np.array(L)
            if sum(x < 1000) != 0:
                results.setdefault(ref, [])
                results[ref] += [list(np.arange(l, l+75)) for l in L]
    f = plt.figure(figsize=(len(results), 5), dpi=300)
    gs = gridspec.GridSpec(len(results), 1)
    gs.update(hspace=0.5)
    pointsd = dict()
    allpoints = []
    for ref in results:
        n = np.array(results[ref])
        n = n[n < 1000]
        np.sort(n)        
        points = Counter(list(n))
        p = []
        for i in range(0, 1000):
            p.append(points[i])
        allpoints += p
        pointsd[ref] = p
    
    answers = getPathogens()
    
    for i, ref in enumerate(pointsd):
        a = plt.subplot(gs[i])
        p = np.array(pointsd[ref])
        if ref.lower() not in answers[str(dataset)]:
            p[p > 10] = 1
        ylim = int(math.ceil(max(p) / 100)) * 100
        a.set_xlim(0, 1001)
        a.set_ylim(0, ylim+1)
        xs = np.arange(0, len(p), 3)
        ys = list(p[xs])
        a.text(0, ylim+(ylim/3), ref.replace("_", " "),
               ha='left', va='top', fontsize=8,
               fontdict={'style':'italic'})
        a.plot([xs, xs], [[0] * len(ys), ys], color='#f0760a')
        a.set_xticks(np.arange(0, 1001, 200))
        a.set_xticklabels(range(0, 1001, 200), fontsize=7)
        a.set_yticks([0, ylim])
        a.set_yticklabels([0, ylim], fontsize=7)
        if i != (len(results) - 1):
            a.xaxis.set_visible(False)
    f.suptitle("Reads Matching Reference Microbe", fontsize=8)
    a.set_xlabel("Position in Microbe Genome", fontsize=7)
    plt.show()
    if save:
        f.savefig("../interface/%s_reads_mapping_%s.png" % (dataset, typ), dpi=150, bbox_inches='tight')
    plt.close()

def showMappingBar(dataset, typ='path', save=False):
    if typ == 'path':
        f = 'references.fasta'
        p = '.pickle'
    elif typ == 'digest':
        f = 'digest.fasta'
        p = '_digest.pickle'
    elif typ == 'env':
        f = 'env.fasta'
        p = '_env.fasta'
    elif typ == 'food':
        f = 'food.fasta'
        p = '_food.fasta' 
    o = open("../data/%s%s" % (dataset, p), "rb")
    r = pickle.load(o)
    o.close()
    bars = []
    labs = []
    for ref, val in r.items():
        bars.append(len(val))
        labs.append(ref)
    f = plt.figure(figsize=(4, len(bars)/4), dpi=300)
    a = f.add_subplot(111)
    a.barh(range(len(bars)), bars, color='#1532a0')
    a.set_yticks(np.arange(0, len(bars)))
    a.set_yticklabels([x.replace("_", " ") for x in labs], fontsize=7)
    a.set_ylim(-0.5, len(bars))
    a.set_xticks(np.arange(0, 4001, 1000))
    a.set_xticklabels(range(0, 4001, 1000))
    a.set_xlabel("Number of Reads", fontsize=8)
    a.set_title("Number of Reads per Reference Microbe", fontsize=8)
    if save:
        f.savefig("../interface/%s_reads_bar_%s.png" % (dataset, typ), dpi=150, bbox_inches='tight')
    plt.show()
    plt.close()


def showMappingTable(dataset, option, typ='path', diet=None):
    if typ == 'path':
        path = "../data/sample_%s_reads.tsv" % dataset
    elif typ == 'food':
        path = "../data/sample_%s_reads_food.tsv" % dataset
    elif typ == 'digest':
        path = '../data/sample_%s_reads_digest.tsv' % dataset
    elif typ == 'env':
        path = '../data/sample_%s_reads_env.tsv' % dataset
    tab = pd.read_csv(path, sep="\t")
    if diet:
        f = pd.read_csv("../data/food_typs.tsv", sep="\t", header=None)
        if diet == 'carnivore':
            f = f[f[1] == 'c']
        elif diet == 'herbivore':
            f = f[f[1] == 'h']
        refnames = list(f[0])
        refnames = [r.capitalize().replace(" ", "_") for r in refnames]
        tab = tab[tab['Reference Microbe'].isin(refnames)]
    if option == 1:
        display(tab)
    elif option == 2:
        webbrowser.open(path)
    elif option == 3:
        return (tab)

def getFoodTable(diet, option=1):
    tab = pd.read_csv("../data/food_viruses.tsv", sep="\t")
    if diet == 'herbivore':
        tab = tab[:18]
    elif diet == 'carnivore':
        tab = tab[18:]
    tab = tab[['virus', 'host']]
    if option == 1:
        display(tab)
    else:
        return (tab)

def plotRefSeq(microbe, save=False):
    nam = microbe
    if nam == "Select option":
        return (None)
    nam = nam.replace(" ", "_")
    nam = nam.replace("_e_", "_E_")
    nam = nam.replace("_a_", "_A_")
    colours = {'A': '#d83b08',
              'C': '#08d80e',
               'T': '#0956e5',
               'G': '#fed400'}
    data = getData()
    lengths = getRefLengths()
    types = dict(zip(data['microbe_name'], data['microbe_type']))
    F = FastaToDict("../data/references.fasta")
    nam2 = nam.lower()
    seq = F[nam]
    f = plt.figure(figsize=(6, 6), dpi=300, edgecolor='black', linewidth=2)
    g = gridspec.GridSpec(5, 4)
    g.update(hspace=0.5, wspace=0.5)

    a = plt.subplot(g[0, 0:2])
    a.text(-0.1, 0.7, "%s" % nam.replace("_", " "), fontsize=12, fontdict={'style': 'italic'})
    a.text(-0.1, 0.4, 'Microbe Type: %s' % types[nam2].replace("es", ""), fontsize=10)
    a.text(-0.1, 0.1, "Microbe Size: %s nucleotides" % splitN(lengths[nam]), fontsize=10)
    a.set_axis_off()
    
    b = plt.subplot(g[1:3, 0:2])
    i = matplotlib.image.imread("../data/images/%s.jpg" % nam)
    b.imshow(i)
    b.set_axis_off()
    
    
    d = plt.subplot(g[1:3, 2:4])
    counts = Counter(list(seq))
    bars = []
    cols = []
    for nuc in ['A', 'C', 'T', 'G']:
        bars.append(counts[nuc] / sum(counts.values()))
        cols.append(colours[nuc])
    d.bar(range(len(bars)), bars, color=cols)
    d.set_xticks(range(len(bars)))
    d.set_xticklabels(['A', 'C', 'T', 'G'], fontsize=4)
    d.set_yticks(np.arange(0, 0.51, 0.1))
    d.set_yticklabels(["%s%%" % x for x in np.arange(0, 51, 10)], fontsize=4)
    d.spines['right'].set_visible(False)
    d.spines['top'].set_visible(False)
    d.yaxis.set_ticks_position('left')
    d.xaxis.set_ticks_position('bottom')
    d.set_title("% of each nucleotide", fontsize=6)
    
    e = plt.subplot(g[3:, 0:])
    tA = re.sub('[GCT]', ' ', seq)
    tA = "\n".join([tA[i:i+100] for i in range(0, 801, 100)])

    tC = re.sub('[GAT]', ' ', seq)
    tC = "\n".join([tC[i:i+100] for i in range(0, 801, 100)])

    tG = re.sub('[ACT]', ' ', seq)
    tG = "\n".join([tG[i:i+100] for i in range(0, 801, 100)])

    tT = re.sub('[GCA]', ' ', seq)
    tT = "\n".join([tT[i:i+100] for i in range(0, 801, 100)])

    e.text(-0.1, 0.9, "First 1000 nucleotides", fontsize=10)
    
    e.text(-0.1, 0.75, ">%s" % nam, fontsize=8,
           fontdict={'fontname': 'monospace'})
    e.text(-0.1, 0, tA, fontsize=8, color=colours['A'],
           fontdict={'fontname': 'monospace', 'weight': 'bold'})
    e.text(-0.1, 0, tC, fontsize=8, color=colours['C'],
           fontdict={'fontname': 'monospace', 'weight': 'bold'})
    e.text(-0.1, 0, tG, fontsize=8, color=colours['G'],
           fontdict={'fontname': 'monospace', 'weight': 'bold'})
    e.text(-0.1, 0, tT, fontsize=8, color=colours['T'],
           fontdict={'fontname': 'monospace', 'weight': 'bold'})
    e.set_axis_off()
    plt.show()
    if save:
        f.savefig("../interface/%s.png" % microbe.replace(" ", "_").lower(), dpi=300, bbox_inches='tight')
    plt.close()

def getMatrixPoints(animal, dataset):
    m2 = pd.read_csv("../data/matrix2.tsv", sep="\t", index_col=0)
    m2['host'] = m2.index.values
    m2 = m2[m2['host'] == animal]
    m2 = m2.drop('host', 1)
    m2.index = ['dataset_%s' % dataset]
    m2.loc[''] = [0] * (len(m2.columns) - 1) + ['#bcbcbc']
    m2 = m2.sort_index()

    return (m2)
    
def plotMatrix(dataset=None):
    matrix = pd.read_csv("../data/matrix_db.txt", sep="\t", index_col=0)
    if dataset:
        dataset_points = getMatrixPoints(getSppAnswers()[str(dataset)], dataset)
        matrix = matrix.append(dataset_points)
    f = plt.figure(figsize=(10, 4), dpi=150)
    a = f.add_subplot(111)
    a.set_ylim(0, len(matrix))

    a.set_xticks(np.arange(0.5, len(matrix.columns)-1, 1))
    a.set_xticklabels([c.replace("_", " ")
    for c in matrix.columns[:-1]], rotation='vertical', fontsize=8)

    a.set_yticks(np.arange(0.5, len(matrix), 1))
    a.set_yticklabels([c.replace("_", " ")
    for c in matrix.index.values], fontsize=8)
    for i, col in enumerate(matrix.columns):
        if col != 'colour':
            for j, row in enumerate(matrix.index.values):
                val = matrix.loc[row, col]
                if val == 1:
                    colour = matrix['colour'].values[j]
                    a.add_patch(matplotlib.patches.Rectangle([i, j], 1, 1,
                                                             color=colour,
                                                             zorder=1))
    a.hlines(np.arange(0, len(matrix)), 0, len(matrix.columns), lw=0.2, zorder=1)
    a.vlines(np.arange(0, len(matrix.columns)), 0, len(matrix), lw=0.2)
    if dataset:
        a.hlines(len(matrix)-2, 0, len(matrix.columns)-1)
        a.hlines(len(matrix)-1, 0, len(matrix.columns)-1)
        a.vlines(np.arange(0, len(matrix.columns)), len(matrix)-2+0.1, len(matrix)-1-0.1,
                 color='white', zorder=10, lw=2)
    a.set_xlim(0, len(matrix.columns)-1)
    a.hlines([5, 15], 0, len(matrix.columns)-1)
    a.text(-5, 17, 'Omnivores', ha='right', color='mediumblue')
    a.text(-5, 12, 'Herbivores', ha='right', color='forestgreen')
    a.text(-5, 3, 'Carnivores', ha='right', color='#fd4a1e')

    plt.show()
    plt.close()

def plotPCA(dataset):
    inf = open("../data/pca_points.pickle", "rb")
    data = pickle.load(inf)
    inf.close()
    answer = getSppAnswers()[str(dataset)]
    colors = ['#fd4a1e',
             '#fd4a1e',
             '#fd4a1e',
             '#fd4a1e',
             '#fd4a1e',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'forestgreen',
             'mediumblue',
             'mediumblue',
             'mediumblue']
    dataset_vals = {'dog': [1.5, 2.2],
                    'pig': [1.5, 1.3],
                    'horse':[ 0.5, 0.5],
                    'cat': [0.5, 2.5],
                    'sheep': [1.4, 0.2],
                    'cow': [1.9, 0.9]}
    f = plt.figure(figsize=(3,3), dpi=300)
    a = f.add_subplot('111')
    for i, key in enumerate(data):
        a.scatter(data[key]['y'], data[key]['x'], facecolor=colors[i], edgecolor='grey', lw=0.2)
        a.text(data[key]['y']+0.015, data[key]['x']+0.07, key.replace("_", " "), fontsize=4, ha='left', va='center')
        a.scatter(dataset_vals[answer][0], dataset_vals[answer][1], facecolor='#fdc71e', edgecolor='black', lw=0.2, marker='*', s=70)
        a.text(dataset_vals[answer][0]+0.015, dataset_vals[answer][1]+0.07,
               "dataset %i" % dataset, fontsize=4, ha='left', va='center')
    a.set_xlim(-0.2, 2.8)
    a.set_ylim(-0.2, 2.8)
    a.set_xticks(np.arange(0, 3, 0.5))
    a.set_yticks(np.arange(0, 3, 0.5))
    a.set_xticklabels(np.arange(0, 3, 0.5), fontsize=6)
    a.set_yticklabels(np.arange(0, 3, 0.5), fontsize=6)
    a.set_ylabel('Similarity', fontsize=6)
    a.set_xlabel('Similarity', fontsize=6)
    plt.show()
    plt.close()


def plotFood(dataset, diet):
    tab1 = getFoodTable(diet, option=3)
    tab2 = showMappingTable(dataset, option=3, typ='food', diet=diet)
    tab1['virus'] = [t.capitalize().replace(" ", "_") for t in tab1['virus']]
    hosts = set(tab1['host'])
    bars = []
    labs = []
    bars2 = []
    labs2 = []
    for host in hosts:
        s1 = tab1[tab1['host'] == host]
        s2 = tab2[tab2['Reference Microbe'].isin(set(s1['virus']))]
        bars.append(sum(s2['Number of Reads in Sample']))
        labs.append(host)
        if sum(s2['Number of Reads in Sample']) != 0:
            bars2.append(sum(s2['Number of Reads in Sample']))
            labs2.append(host)
    f = plt.figure(figsize=(6, 3), dpi=150)
    a = f.add_subplot(121)
    cmap = matplotlib.cm.Paired(np.linspace(0, 1, 6))
    a.bar(range(len(bars)), bars, color=cmap)
    a.set_xticks(np.arange(len(bars)))
    a.set_xticklabels(labs, rotation='vertical')
    a.set_xlabel('Food')
    a.set_ylabel("Number of Reads")
    b = f.add_subplot(122)
    b.pie(bars2, colors=cmap, labels=labs2, autopct='%.2f%%')
    plt.show()
    plt.close()