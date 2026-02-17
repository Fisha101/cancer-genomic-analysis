import os
import csv #for report
import matplotlib.pyplot as plt 

genes = []                       #for graph
gc_values = []

report = open('report.csv', 'w', newline='')
writer = csv.writer(report)
writer.writerow(['Gene', 'Length_bp', 'GC%', 'Longest_ORF_bp'])

files = sorted(os.listdir('data'))   # <- optional but nice

for i in files:
    if i.endswith('.fasta'):
        print("opening:", i)
        f = open("data/" + i)
        sequence = ''
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line
        sequence = sequence.upper()

        no_c = sequence.count('C')
        no_g = sequence.count('G')
        gc_percent = (no_c + no_g)*100/len(sequence)

        gene = i.replace('.fasta','')         
        genes.append(gene)
        gc_values.append(round(gc_percent, 2))

        stops = ['TAA', 'TAG', 'TGA']
        longest_orf = 0
        for frame in [0,1,2]: 
            for j in range(frame, len(sequence)- 2, 3):
                codon = sequence[j:j+3]
                if codon == 'ATG':
                    for k in range(j+3,len(sequence)-2,3):
                        if sequence[k:k+3] in stops:
                            current_length = (k+3) - j
                            if current_length > longest_orf:
                                longest_orf = current_length
                            break                

        print(gene)
        print('Length:', len(sequence), 'bp')
        print('GC%:', round(gc_percent, 2), '%')
        print('Longest ORF:', longest_orf, 'bp')
        print("-"*30)

        writer.writerow([gene, len(sequence), round(gc_percent, 2), longest_orf])
        f.close()

report.close()

# plot GC% 
plt.bar(genes, gc_values)
plt.ylabel("GC %")
plt.title("GC content by gene")
plt.savefig("plots/gc_content.png")

