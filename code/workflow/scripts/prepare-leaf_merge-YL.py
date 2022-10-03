
from subproc import run

subsets = run("ls ../run/full/subset_*_refined")[0].split()

N = 100

fout = open("run_merge.sh",'w')
for i in range(0,len(subsets)/N+1):
    libs = subsets[i*N:(i+1)*N]

    fout.write("sbatch -t 12:00:00 --partition=broadwl --mem 4G -o log --wrap=\'python2 merge_leafcutter_clusters.py ../run/full/merged_%d 3 %s\'\n"%(i," ".join(libs)))

fout.close()

print "python2 merge_leafcutter_clusters.py ../run/full/full_merged_refined 3 %s"%(" ".join(["../run/full/merged_%d"%x for x in range(0,len(subsets)/N+1)]))

subsets = run("ls ../run/full/subset_*_noisy_introns")[0].split()

fout = open("run_merge_noisy.sh",'w')
for i in range(0,len(subsets)/N+1):
    libs = subsets[i*N:(i+1)*N]
    fout.write("sbatch -t 12:00:00 --partition=broadwl --mem 4G -o log --wrap=\'python2 merge_leafcutter_clusters.py run/full/mergednoisy_%d 3 %s\'\n"%(i," ".join(libs)))

fout.close()

#print "python2 merge_leafcutter_clusters.py ../run/full/full_merged_refined %s"%(" ".join(["../run/full/merged_%d"%x for x in range(0,len(subsets)/N+1)]))

print "python2 merge_leafcutter_clusters.py ../run/full/full_mergednoisy_refined 3 %s"%(" ".join(["../run/full/mergednoisy_%d"%x for x in range(0,len(subsets)/N+1)]))

subsets = run("ls ../run/full/subset_*_summaryusage")[0].split()


N = 50
fout = open("run_merge_summary.sh",'w')
for i in range(0,len(subsets)/N+1):
    libs = subsets[i*N:(i+1)*N]
    fout.write("sbatch -t 12:00:00 --partition=broadwl --mem 6G -o log --wrap=\'python2 parse_summaryusage.py run/full/merge_%d_summary %s\'\n"%(i," ".join(libs)))

fout.close()
print "python2 parse_summaryusage.py run/full/full_summary %s"%(" ".join(["../run/full/merge_%d_summary"%x for x in range(0,len(subsets)/N+1)]))



N = 50
subsets = run("ls run/full/subset_*_run/subset_*_perind.constcounts.gz")[0].split()

fout = open("run_merge_cryptic.sh",'w')
for i in range(0,len(subsets)/N+1):
    libs = subsets[i*N:(i+1)*N]
    fout.write("sbatch -t 2:00:00 --partition=broadwl --mem 16G -o log --wrap=\'python2 get_cryptic_counts.py run/full/cryptic_merged_%i.pck %s\'\n"%(i," ".join(libs)))

fout.close()



fout = open("run_merge_table.sh",'w')
for i in range(0,len(subsets)/N+1):
    libs = subsets[i*N:(i+1)*N]
    fout.write("sbatch -t 36:00:00 --partition=broadwl --mem 16G -o log --wrap=\'python2 merge_table.py run/full/merged_%d_perind.constcounts.gz %s\'\n"%(i," ".join(libs)))

fout.close()


print "python2 merge_table.py run/full/All_perind.constcounts.gz run/full/merged_*_perind.constcounts.gz"

#print "python2 parse_summaryusage.py run/full/full_summary %s"%(" ".join(["../run/full/merge_%d_summary"%x for x in range(0,len(subsets)/N+1)]))
