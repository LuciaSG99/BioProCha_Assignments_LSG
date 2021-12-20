require 'bio'
s_pombe_fa,a_thaliana_fa = ARGV
### take protein X in Species A, and BLAST it against all proteins in Species B.  The top hit in Species B is then BLASTed against all proteins in Species A.
##If it’s top hit is the Protein X, then those two proteins are considered to be Orthologue candidates. 
##Using BioRuby to blast and parse the blast reports, find the orthologue pairs between species Arabidopsis and S. pombe.
##I have uploaded their complete proteomes to the Moodle for you. the existing BioRuby objects are sufficient

#CAREFULL!! The fa file of S.Pombe contains aminoacid sequences and the A.Thaliana contains nucleotide sequences

#Things I have to do:

# 1) Create the two databases
# 2) Parse one database to use each sequence as the first query (protein X for example of S.Pombe) and performance the blast against all proteins in A.Thaliana.
#As S.Pombe sequence are aa, in this case I must performance tblastn (TBLASTN search translated nucleotide databases using a protein query).
#Of the output, choose and store the one that has better parameters (best hit, best evalue and other¿?)
# 3) Them do a blastx (BLASTX search protein databases using a translated nucleotide query) of that hit (protein Y) against the database of S.Pombe. Again store the best hit 
# 4) If this best hit == protein X --> Prot Y and X are orthologs. Store this pair in a dictionary? array?
# 5) Do a report with this pairs

#For 1,2,3 and 4 design just one function?? 

######### 1.CREATE THE DATABASES:

def create_database(name_file,type_db,name_db)
 system("makeblastdb -in #{name_file} -dbtype #{type_db} -out #{name_db}")#using the function "system" to run a command in the shell:
end

create_database(s_pombe_fa,'prot','s_pombe_db')
create_database(a_thaliana_fa,'nucl','a_thaliana_db')

#create the S.Pombe BLAST database. Because it has sequences of aminoacids, we must use the method blastx with the nucleotides queries of A.thaliana:
s_pombe_db = Bio::Blast.local('blastx','s_pombe_db') 

#create the A.Thaliana BLAST database. Because it has sequences of nucleotides, we must use the method tblastn with the aminoacides queries of S.Pombe:
a_thaliana_db = Bio::Blast.local('tblastn','a_thaliana_db') 

######### 2. PERFORMANCE THE "RECIPROCAL BEST BLAST” TO FIND PUTATIVE ORTHOLOGUES:
### FILTER PARAMETERS:
# maximum E-value threshold of 1×10 − 6
# coverage of at least 50%

##REFERENCES:
# 1. Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, Bioinformatics, Volume 24, Issue 3, 1 February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585
# 2. Ward, N., & Moreno-Hagelsieb, G. (2014). Quickly finding orthologs as reciprocal best hits with BLAT, LAST, and UBLAST: how much do we miss?. PloS one, 9(7), e101850. https://doi.org/10.1371/journal.pone.0101850

s_pombe_sequences = Bio::FlatFile.open(s_pombe_fa) #flat file to contain the protein sequences to use each as a query
a_thaliana_sequences = Bio::FlatFile.open(a_thaliana_fa) #flat file to contain the protein sequences to use each as a query

#s_pombe_sequences.each{|query_pombe| reciprocal_best_blast(query_pombe,a_thaliana_db,s_pombe_db,1e-6,0.5)}
query_1 = s_pombe_sequences.entries[23]
puts "query inspect: #{query_1.entry_id}"
putative_orthologues = Hash.new #Hash that will contain the orthologues discovered


def reciprocal_best_blast(query1, database1, database2,evalue_thresh,coverage_thres)
  
  ##1. Perform the first blast search:
  results_1 = database1.query(query1) 
  best_hit_1 = results_1.hits[0] #blast orders the results be the evalue, starting from the lowest value, so the first hit is the one with the best evalue
  
  coverage_hit = (best_hit_1.query_end.to_f - best_hit_1.query_start.to_f)/(best_hit_1.query_len.to_f) # calculate the coverage proportion of the alignment
  
  if best_hit_1.evalue.to_f <= evalue_thresh.to_f && coverage_hit.to_f >= coverage_thres.to_f #if the first hit satisfy the filter that was established:
     
     
        
       puts "\ninspect hit: #{best_hit_1.definition.split("|")}"
       ## 2. Perform the second blast search --> Here best_hit_1 is used as a query against the second database:
       results_2 = database2.query(best_hit_1.target_seq)
       best_hit_2 = results_2.hits[0]
       #puts "\nquery1_id: #{query1.definition}"
       #puts "\nhit_second_blast_id: #{best_hit_2.target_id}"
       #puts "\nsequence query 1: #{query1.seq}"
       #puts "\nsequence second hit: #{best_hit_2.target_seq}"
    
     
  else
   puts "\nevalue: #{best_hit_1.evalue} and coverage: #{coverage_hit}"
  end
      
end

reciprocal_best_blast(query_1,a_thaliana_db,s_pombe_db,10e-6,0.5)
#s_pombe_sequences.each {|entry|  list_seq_pombe_prueba << entry} #puts "ID sequence: #{entry.entry_id}"
    #puts "sequence: #{entry.seq}"}
#puts list_seq_pombe_prueba[0..2].each {|query| puts query.seq}  # para ver las queries



##### To see the output:
#results.each {|hit| puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} at "
#  puts "coverage = #{hit.cov}"
#  puts "#{hit.lap_at}"   # this tells you start and end of both the query and the hit sequences
#  hit.each do |hsp|
#    puts hsp.qseq  # this is the gapped Alignment as text of the query
#    puts hsp.hseq
#end}

### Store the best hit:
#hits = Array.new
#results.each {|hit| hits << hit} #store the hits into an Array
#
#best_hit = results[0] #BLAST always order the hits, starting from the lowest evalue

#puts "start of sequence align: #{best_hit_1.query_start.to_i}"
  #puts "end of sequence align: #{best_hit_1.query_end.to_i}"
  #puts "length of query sequence: #{best_hit_1.query_len.to_i}"

