from itertools import islice
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='''Spectral library feature extractor''')

parser.add_argument('infile', metavar='-ip', type=str, nargs='+', help='OpenSWATH compatible spectral library in tsv format')

args = parser.parse_args()

#infile = "OSCC_Crude_Chewers_Smokers_Multi-Consesus_Spec_Lib.tsv"

def get_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            pep = split_i.index("PeptideSequence")
            mod_pep = split_i.index("ModifiedPeptideSequence")
            pro = split_i.index("ProteinId")
            mz = split_i.index("PrecursorMz")
            z = split_i.index("PrecursorCharge")
            ms2_mz = split_i.index("ProductMz")
            ms2_z = split_i.index("ProductCharge")
            frag = split_i.index("FragmentIonType")
            return pep, mod_pep, pro, mz, z, ms2_mz, ms2_z, frag

def extract_dicts(indict):
    sum_total = 0
    for k, values in indict.items():
        peps = {v: v for v in values}
        sum_total += len(peps)

    try:
        dicts_keys = sorted({int(k) : v for k, v in indict.items()})
    except:
        dicts_keys = sorted(indict.keys())#, key=lambda item: int(item[0]))
        
    output = []
    for key in dicts_keys:
        if str(key) in indict: 
            peps = {v: v for v in indict[str(key)]}
            output.append([str(key), str(len(peps)), str(len(peps)/sum_total*100)])
    return output

def count_dicts_keys(indict):
    count_keys = {}
    for key, values in indict.items():
        if len(values) not in count_keys:
            count_keys[len(values)] = [key]
        else:   
            count_keys[len(values)].append(key)
    try:
        dicts_keys = sorted({int(k) : v for k, v in count_keys.items()})
    except:
        dicts_keys = sorted(count_keys.keys())#, key=lambda item: int(item[0]))

    sum_total = 0
    for k, values in count_keys.items():
        sum_total += len(values)

    output = []  
    for key in dicts_keys:
        if key in count_keys:
            output.append([str(key), str(len(count_keys[key])), str(len(count_keys[key])/sum_total*100)])    
    return output


def plot_bar(infile, outfile_suffix, output_list, columns_list):   
    plot_pdf = infile.rstrip(".tsv") + "_" + str(outfile_suffix) + "_barplot.pdf"
    plot_png = infile.rstrip(".tsv") + "_" + str(outfile_suffix) + "_barplot.png"
    
    df = pd.DataFrame(output_list, columns = columns_list, dtype = int)
    df[[columns_list[0], columns_list[1]]] = df[[columns_list[0], columns_list[1]]].astype(int)
    #sns.set_theme(style="white")
    sns.barplot(x=columns_list[0], y=columns_list[1], data = df)
    sns.despine()
    plt.xlabel(columns_list[0],fontsize=10)
    plt.ylabel(columns_list[1],fontsize=10)
    plt.tick_params(labelsize=5)
    plt.savefig(plot_pdf, dpi=300, orientation="landscape")
    plt.savefig(plot_png, dpi=300, orientation="landscape")

def plot_cat_bar(infile, outfile_suffix, output_list, columns_list):  
    plot_pdf = infile.rstrip(".tsv") + "_" + str(outfile_suffix) + "_barplot.pdf"
    plot_png = infile.rstrip(".tsv") + "_" + str(outfile_suffix) + "_barplot.png"
    df = pd.DataFrame(output_list, columns = columns_list, dtype = int)
    
    df[[columns_list[0]]] = df[[columns_list[0]]].astype('category')
    df[[columns_list[1]]] = df[[columns_list[1]]].astype(int)
    #sns.set_theme(style="white")
    sns.catplot(x=columns_list[0], y=columns_list[1], kind="bar", data = df)
    sns.despine()
    plt.xlabel(columns_list[0],fontsize=10)
    plt.ylabel(columns_list[1],fontsize=10)
    plt.tick_params(labelsize=5)
    plt.savefig(plot_pdf, dpi=300, orientation="landscape")
    plt.savefig(plot_png, dpi=300, orientation="landscape")


def extract_features(infile):        
    a = get_header_idx(infile)
    pep_lens = {}
    pep_z = {}
    frag_z = {}
    frags_sum = {}
    frag_ions = {}
    pep_frags = {}
    pro_peps = {}
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            pep_len = str(len(split_i[a[0]]))
            if pep_len not in pep_lens:
                pep_lens[pep_len] = [split_i[a[1]]]
            else:
                pep_lens[pep_len].append(split_i[a[1]])
            if split_i[a[4]] not in pep_z:
                pep_z[split_i[a[4]]] = [split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]]
            else:
                pep_z[split_i[a[4]]].append(split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]])

            if split_i[a[6]] not in frag_z:
                frag_z[split_i[a[6]]] = [split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]] + "@" + split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]]
            else:
                frag_z[split_i[a[6]]].append(split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]] + "@" + split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]])

            if split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]] not in frags_sum:
                frags_sum[split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]] = [split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]]]
            else:
                frags_sum[split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]].append(split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]])

            if split_i[a[-1]] not in frag_ions:
                frag_ions[split_i[a[-1]]] = [split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]] + "@" + split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]]
            else:
                frag_ions[split_i[a[-1]]].append(split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]] + "@" + split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]])

            if split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]] not in pep_frags:
                pep_frags[split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]] = [split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]]]
            else:
                pep_frags[split_i[a[1]] + "@" + split_i[a[3]] + "@" + split_i[a[4]]].append(split_i[a[-1]] + "@" + split_i[a[6]] + "@" + split_i[a[5]])

            if split_i[a[2]] not in pro_peps:
                pro_peps[split_i[a[2]]] = [split_i[a[1]]]
            else:
                pro_peps[split_i[a[2]]].append(split_i[a[1]])
                

    pep_lens_output = extract_dicts(pep_lens) 
    pep_z_output = extract_dicts(pep_z)
    frag_z_output = extract_dicts(frag_z)
    frag_ions_output = extract_dicts(frag_ions)
    pro_peps_output = extract_dicts(pro_peps)
    pep_frags_output = count_dicts_keys(pep_frags)

    pep_lens_outfile = "{0}SpecLib_peptide_length_summary.txt".format(infile.rstrip(".tsv"))
    with open(pep_lens_outfile, "w") as outf:
        outf.write("Peptide length\tNo. of Peptides\tPercentage(%)\n")
        outf.writelines("\t".join(i) + '\n' for i in pep_lens_output)

    pep_z_outfile = "{0}SpecLib_peptide_chargestate_summary.txt".format(infile.rstrip(".tsv"))
    with open(pep_z_outfile, "w") as outf:
        outf.write("Precursor charge\tNo. of Peptides\tPercentage(%)\n")
        outf.writelines("\t".join(i) + '\n' for i in pep_z_output)

    frag_z_outfile = "{0}SpecLib_fragment_chargestate_summary.txt".format(infile.rstrip(".tsv"))
    with open(frag_z_outfile, "w") as outf:
        outf.write("Fragment charge\tNo. of Fragments\tPercentage(%)\n")
        outf.writelines("\t".join(i) + '\n' for i in frag_z_output)

    frag_ions_outfile = "{0}SpecLib_fragment_ions_summary.txt".format(infile.rstrip(".tsv"))
    with open(frag_ions_outfile, "w") as outf:
        outf.write("Fragment Ion Type\tNo. of Fragments\tPercentage(%)\n")
        outf.writelines("\t".join(i) + '\n' for i in frag_ions_output)

    pep_frags_outfile = "{0}SpecLib_frag_ions_per_precursor_summary.txt".format(infile.rstrip(".tsv"))
    with open(pep_frags_outfile, "w") as outf:
        outf.write("Fragment Ions\tNo. of Precursors\tPercentage(%)\n")
        outf.writelines("\t".join(i) + '\n' for i in pep_frags_output)

    pro_peps_outfile = "{0}SpecLib_peptides_per_proteins_summary.txt".format(infile.rstrip(".tsv"))
    with open(pro_peps_outfile, "w") as outf:
        outf.write("Proteins\tNo. of Peptides\tPercentage(%)\n")
        outf.writelines("\t".join(i) + '\n' for i in pro_peps_output)

    plot_bar(infile, "peptide_length", pep_lens_output, ["Peptide length", "No. of Peptides", "Percentage(%)"])
    plot_bar(infile, "peptide_chargestate", pep_z_output, ["Precursor charge", "No. of Peptides", "Percentage(%)"])
    plot_bar(infile, "fragment_chargestate", frag_z_output, ["Fragment charge", "No. of Fragments", "Percentage(%)"])
    plot_cat_bar(infile, "fragment_ions", frag_ions_output, ["Fragment Ion Type", "No. of Fragments", "Percentage(%)"])
    plot_bar(infile, "frag_ions_per_precursor", pep_frags_output, ["Fragment Ions", "No. of Precursors", "Percentage(%)"])

if __name__== "__main__":
    extract_features(args.infile[0])
