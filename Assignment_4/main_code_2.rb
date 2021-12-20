require 'bio'
s_pombe_fa,a_thaliana_fa = ARGV #store the .fa files in variables

#CAREFULL!! The .fa file of S.Pombe contains aminoacid sequences and the A.Thaliana contains nucleotide sequences

################################################## 1.CREATE THE DATABASES ##############################################################################

#Function that creates a blast database by a .fa file, type of database we want and name of the database as arguments:
def create_database(name_file,type_db,name_db)
 system("makeblastdb -in #{name_file} -dbtype #{type_db} -out #{name_db}")#using the function "system" to run a command in the shell:
end

create_database(s_pombe_fa,'prot','s_pombe_db') #apply the function using the .fa file of S.Pombe, which has sequences of aminoacids
create_database(a_thaliana_fa,'nucl','a_thaliana_db')#apply the function using the .fa file of A.Thaliana, which has sequences of nucleotides

#create the S.Pombe BLAST database. Because it has sequences of aminoacids, we must use the method blastx with the nucleotides queries of A.thaliana:
s_pombe_db = Bio::Blast.local('blastx','s_pombe_db') 

#create the A.Thaliana BLAST database. Because it has sequences of nucleotides, we must use the method tblastn with the aminoacides queries of S.Pombe:
a_thaliana_db = Bio::Blast.local('tblastn','a_thaliana_db') 

################################################ 2. PERFORMANCE THE "RECIPROCAL BEST BLAST” TO FIND PUTATIVE ORTHOLOGUES ##################################
### FILTER PARAMETERS:
# maximum E-value threshold of 1×10 − 6
# coverage of at least 50%

## REFERENCES:
# 1. Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, Bioinformatics, Volume 24, Issue 3, 1 February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585
# 2. Ward, N., & Moreno-Hagelsieb, G. (2014). Quickly finding orthologs as reciprocal best hits with BLAT, LAST, and UBLAST: how much do we miss?. PloS one, 9(7), e101850. https://doi.org/10.1371/journal.pone.0101850

s_pombe_sequences = Bio::FlatFile.open(s_pombe_fa) #flat file that contains the protein sequences to use each as a query
a_thaliana_sequences = Bio::FlatFile.open(a_thaliana_fa) #flat file that contains the nucleotides sequences to use each as a query


putative_orthologues = [] #Array that will contain the orthologues discovered

# FUNCTION TO OBTAIN THE BEST HIT OF A BLAST SEARCH: As arguments it needs a sequence query, a database, a evalue and a coverage threshold to use as a filter
def obtain_best_hit(query, database,evalue_thresh,coverage_thres)
  
  ##1. Perform the blast search:
  results = database.query(query) #store the results of the blast alignment
  best_hit = results.hits[0] #blast orders the results by the evalue, starting from the lowest value, so the first hit is the one with the best evalue
  
  if best_hit == nil #if there are not any hit:
   
    return "NA"
  
  else
    coverage_hit = (best_hit.query_end.to_f - best_hit.query_start.to_f)/(best_hit.query_len.to_f) # calculate the coverage proportion of the alignment
    
    if best_hit.evalue.to_f <= evalue_thresh.to_f && coverage_hit.to_f >= coverage_thres.to_f #if the first hit satisfies the filter that was established:
    
      return best_hit #return that hit as a results of the function 
    
    else #if not:
    
     return "NA" #return a NA (string) as a value
    end
  end  
end

### APPLYING THE FUNCTION TWICE TO SEARCH FOR PUTATIVE ORTHOLOGUES:
s_pombe_sequences.entries.each {|entry_pombe| #for each sequence of the S.Pombe file:
 
 ##### 1. Perform the first blast search of the sequence of Pombe (aa) against the A.Thaliana Database (nuclt)
  best_hit_1 = obtain_best_hit(entry_pombe,a_thaliana_db,10e-6,0.5) #apply the function with the query, the A.Thaliana database, with a evalue and coverage threshold:
  
  if best_hit_1 == "NA" #if there aren't hits skip to the next query
    next
  else
   
    #Obtain the best hit of the A.Thaliana as a query:
    a_thaliana_sequences.rewind() #reset to search in the flat file from the beginning
    a_thaliana_sequences.entries.each{|entry_thaliana| #for each sequence of the A.Thaliana file
    
    if best_hit_1.definition.split("|")[0].rstrip == entry_thaliana.entry_id # to obtain that hit as a query in from the flat file that contains the sequences:
    
      ##### 2. Perform the second blast search of the sequence of Thaliana (nuclt) against the S.Pombe Database (aa)
      best_hit_2 = obtain_best_hit(entry_thaliana,s_pombe_db,10e-6,0.5) #apply the function again with the new query
    
      unless best_hit_2 == "NA" #unlesss there is no hits:
        if best_hit_2.definition.split("|")[0].rstrip == entry_pombe.entry_id # if they are putative orthologues (a.k.a. this hit matches with the first query:
          puts "FOUND PUTATIVE ORTHOLOGUES!\n"
          array_pair = [entry_thaliana.entry_id,entry_pombe.entry_id] #store in an array the sequences that are orthologues. 
          putative_orthologues << array_pair
         
        else
          next
    
        end
      end  
    end}
         
  end
}
   
################################################## 3. CREATE THE REPORT ########################################################################################

report = File.open("Report_putatives_orthologues.txt", "w")

report << 'REPORT PUTATIVE ORTHOLOGUES BY APPLYING THE METHOD --> RECIPROCAL BEST BLAST'
report << "\n\nThere were found #{putative_orthologues.length} pairs of putative orthologues doing a Reciprocal best blast using as a filter an evalue <= 1e-6 and a minimum coverage of 50%"
report << "\n\nReferences used to establish the parameters:\n"
report << "\n1. Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, Bioinformatics, Volume 24, Issue 3, 1 February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585"
report << "\n2. Ward, N., & Moreno-Hagelsieb, G. (2014). Quickly finding orthologs as reciprocal best hits with BLAT, LAST, and UBLAST: how much do we miss?. PloS one, 9(7), e101850. https://doi.org/10.1371/journal.pone.0101850\n"
report << "\n\n\nID_A.THALIANA\t\tID_S.POMBE\n"
putative_orthologues.each{|pair| report << "\n#{pair[0]}\t\t#{pair[1]}"}

report.close