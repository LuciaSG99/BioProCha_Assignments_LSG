### ASSIGNMENT 3 : MAIN CODE
require 'rest-client'  
require 'bio'
list_genes = ARGV[0] #obtain the list of genes as a argument from the shell

######## 1. RETRIEVE THE SEQUENCES OF THE GENES: 

#Function that saves the genes of the list in an array:
def read_file_genes(file)
  genes_list = Array.new #array that will contain the genes names
  lines = IO.readlines(file)
  lines.each {|gene| genes_list << gene.chomp} #chomp to remove the "\n" of each line
  return genes_list #return the array with the genes
end

genes = read_file_genes(list_genes)

#Function that creates a EMBL object by a gene name
def create_obj_embl(gene)
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}"  
  response = RestClient.get(address)
  record = response.body #store the result of the search
  entry = Bio::EMBL.new(record) #create a new EMBL object
  return entry 
end

#Creation of a Hash which will contain the genes of the list as a EMBL objects from BioRuby:
hash_embl_objs = Hash.new #hash that will contain the EMBL objects, as the keys the gene name 
genes.each {|gene| hash_embl_objs[gene.upcase] = create_obj_embl(gene)} #apply the function


########## 2 AND 3. SEARCH FOR THE CTTCTT SEQUENCE IN EACH EXON AND CREATE NEW SEQUENCE FEATURE:

#FUNCTION TO OBTAIN THE POSITION OF THE SEQUENCE IN THE GENE OR CHROMOSOME:

def obtain_position_cttctt(key,embl_obj,feature,exon,chr_start) #as arguments: the gene name, the EMBL object, the feature, exon == 1 argument if we want the position in the gene
 #chr_starts is the first position of the chromosome in case we want the sequence position in the whole chromosome coordenates
  re = Regexp.new(key)
  position_matches = Array.new() #Array that will contain the coordenates of all the matches for one gene
  
  if feature.feature == "exon" #if the feature is a exon:
     
      if feature.qualifiers[0].value.match(re)#if it is a exon that belongs to the gene of the list:
      
        if feature.position.match(/complement\(\d+\.\.\d+\)/) == nil  and feature.position.match(/\d+\.\.\d+/) == nil #ŧo avoid repeated genes or some other elements 
            return
        end
        
        localization_obj = feature.locations.locations[0] #obtain the information about the position of the exon
        start = localization_obj.from.to_s #obtain the start position of the exon                
        seq_exon = embl_obj.seq.splice(feature.position) #store the sequence of the exon
        
        if localization_obj.strand == 1 #if the exon is in the forward strand:            
                      
            matches = seq_exon.seq.to_enum(:scan,/(?=(cttctt))/).map{Regexp.last_match} #to get all the matches in overlapping situations like this: CTTCTTCTT 
            
            unless matches.empty? #unless no matches were found:      
               
              matches.each{|match|
              
            ###Obtain the coordenates of the match in the whole gene sequence:
            
              if exon == 1 #in case we want the coordenates in the gene sequence:
               
               #to obtain the coordenates of the match in the whole sequence, we have to sum the initial position (nucleotide) of the exon plus the initial
               #position of the match:
                start_pos_seq = match.begin(0).to_i  + localization_obj.from.to_i #initial position of the match in the sequence (GENE)
                
              else #if we want the chromosomes coordenates of the sequence, we have to sum the start position of the match, of the exon and of the chromosome:
               
                start_pos_seq = match.begin(0).to_i  + start.to_i + chr_start.to_i - 1.to_i #initial position of the match in the whole chromosome sequence 
              end
                                           
              end_pos_seq = start_pos_seq.to_i + 5.to_i #the final position of the match is calculated by adding 5 to the start position (because the repeat has 6 nucleotides)
              
              coordenates = "#{start_pos_seq.to_s}..#{end_pos_seq.to_s}" 
                           
              position_matches << coordenates} #insert the coordenates into the array
              
            end
        
        end  
            
        if localization_obj.strand == -1 #if the exon is in the complementary reverse strand
                    
            #Because the exon is in the complementary reverse strand, we can obtain the position of the match CTTCTT by searching for the sequence "AAGAAG" in
            #the forward strand. 
           matches_2 = seq_exon.seq.to_enum(:scan,/(?=(aagaag))/).map{Regexp.last_match}  
                      
           unless matches_2.empty? #unless no matches were found:      
                      
             matches_2.each{|match|
             if exon == 1
                start_pos_seq = match.begin(0).to_i  + localization_obj.from.to_i #initial position of the match in the sequence (GENE)
              
             else
              
                start_pos_seq = match.begin(0).to_i  + start.to_i + chr_start.to_i - 1.to_i #initial position of the match in the whole chromosome sequence 
                
             end
                                           
              end_pos_seq = start_pos_seq.to_i + 5.to_i #the final position of the match is calculated by adding 5 to the start position (because the repeat has 6 nucleotides)
              coordenates = "complement(#{start_pos_seq.to_s}..#{end_pos_seq.to_s})"
              
              position_matches << coordenates}
           end
        end
      end
  end
  
  return position_matches #return the array of the coordenates of each match in the exon         
end

#FUNCTION TO CREATE NEW FEATURES : 
def new_feature(name,coordenates,strand,id_gene_exon) #name = name of the feature, coordenates of the match, id = gene and exon where the match was found
  feature = Bio::Feature.new(feature = name,position = coordenates) #make a new feature
  feature.append(Bio::Feature::Qualifier.new('strand', strand)) #indicates in what strand it is 
  feature.append(Bio::Feature::Qualifier.new('id', id_gene_exon)) #Gene_Name_Exon_number as an id
  return feature  #return the feature         
end

########Store the exons sequences for each gene and search for the repetitions in the exons:
hash_bio_seq_obj = Hash.new #to add the new features we have to use Bio::Sequence objects instead of Bio:EMBL objects. This hash will contain those objects 
hash_embl_objs.each{|gene_name,emblobj| #for each key (Gene name, for example = "AT4g27030") and EMBL objects:
  bio_seq_obj = emblobj.to_biosequence #convert the EMBL object into a Biosequence object to introduce the features
  emblobj.features.each {|ft| #for each exon:
   #transforms the EMBL object into a Bio sequence object
  list_coordenates_gene = obtain_position_cttctt(gene_name,emblobj,ft,1.to_i,"nothing") #apply the function to obtain the coordenates
  
  list_coordenates_gene.each {|coordenates| #for each match coordenates
   
   if ft.locations.locations[0].strand == -1 #if it is in the reverse complementary strand:
    
    feature_match = new_feature("CTTCTT_match",coordenates,'-',"#{gene_name}_#{ft.qualifiers[0].value.match(/exon\d/)[0]}")#apply the function
   end
   
   if ft.locations.locations[0].strand == +1 #if it is in the forward strand:
    
    feature_match = new_feature("CTTCTT_match",coordenates,'+',"#{gene_name}_#{ft.qualifiers[0].value.match(/exon\d/)[0]}")
       
   end
       bio_seq_obj.features << feature_match  #insert the new feature in the Biosequence object
      
  }}
   
  hash_bio_seq_obj[gene_name] = bio_seq_obj} #insert the Biosequence object in the hash


###### Search for genes that do not contain the sequence in their exons:
list_genes_no_repeat = Array.new #Array that will contain those genes

 hash_bio_seq_obj.each{|gene_report,bio_obj| #for each gene and its Biosequence object:
  #if that gene does not have features called CTTCTT_match, it means that do not have that sequence
  #in its exons because that type of feature was only created when there was a match
  if bio_obj.features.any?{|row| row.feature == "CTTCTT_match"} == false 
   list_genes_no_repeat << gene_report #store the name of the gene in the array
  end}

#### Store only the genes that have the sequence CTTCTT at least in one: 
hash_bio_obj_rep = Hash.new #hash that will contain only those genes with the repeated sequence 
hash_bio_seq_obj.each{|gene_repeat,obj|
 unless list_genes_no_repeat.include? gene_repeat 
  hash_bio_obj_rep[gene_repeat] = obj 
 end}

  
########### 4A. CREATING A GFF3-FORMATTED FILE OF THE FEATURES:

gff3_file = File.open("new_features.gff", 'w') 
gff3_file << "##gff-version 3\n\n"
hash_bio_obj_rep.each{|g_n,value|
 value.features.each {|featu|
  if featu.feature == "CTTCTT_match" #if it is one of the features we create
    if featu.position.match(/complement/) 
     coordenates = featu.position.split("(") #in case that the match is in the complement reverse strand
     coordenates = coordenates[1].chomp(")")
     start_posit,end_posit = coordenates.split("..")
    
    else 
     start_posit,end_posit = featu.position.split("..")
    end
    gff3_file << "#{g_n}\tBioRuby\trepeated_sequence\t#{start_posit}\t#{end_posit}\t.\t#{featu.assoc["strand"]}\t.\t#{featu.assoc["id"]}\n\n"
  end}
  
  }
gff3_file.close

############. 4B. REPORT OF GENES WITHOUT SEQUENCE CTTCTT:
## list_genes_no_repeat was obtained previously 
report_file = File.open("Report_genes_no_repeat.txt","w")
report_file << "REPORT OF GENES THAT DO NOT HAVE THE SEQUENCE 'CTTCTT' IN THEIR EXONS:\n\n"
report_file << "Number of total genes from the list: #{genes.length}\n"
report_file << "Number of genes without the sequence CTTCTT: #{list_genes_no_repeat.length}\n\n"
list_genes_no_repeat.each {|gene| report_file << gene
 report_file << "\n"}


############ 5.  GFF3 file WITH THE FULL CHROMOSOME COORDINATES:
#Obtain chromosomes start positions (for then apply the function obtain_position_cttctt):

hash_bio_seq_obj = Hash.new #to add the new features we have to use Bio::Sequence objects with the chromosome coordenates
hash_bio_obj_rep.each{|nm,obj| #for each key (Gene name, for example = "AT4g27030") and Biosequence object
  start_chr = obj.primary_accession.split(":")[3] #the start position of the chr
  
  obj.features.each {|ft2| 
  list_coordenates_gene = obtain_position_cttctt(nm,obj,ft2,0,start_chr) #because we want the chromosomes coordenates, we should put any number except 1 and give the chromosome start position
  list_coordenates_gene.each {|coordenates| 
   if ft2.locations.locations[0].strand == -1 #if it is the reverse complementary strand
    feature_match = new_feature("CTTCTT_match_chr",coordenates,'-',"#{nm}_#{ft2.qualifiers[0].value.match(/exon\d/)[0]}")
   end
   
   if ft2.locations.locations[0].strand == +1 #if it is the forward strand
    feature_match = new_feature("CTTCTT_match_chr",coordenates,'+',"#{nm}_#{ft2.qualifiers[0].value.match(/exon\d/)[0]}")
   end
   obj.features << feature_match #insert the new feature 
   
   }}
   
  hash_bio_seq_obj[nm] = obj}


#### WRITE THE GFF3 FILE:
gff3_file = File.open("new_features_chr.gff", 'w')
gff3_file << "##gff-version 3\n\n"
hash_bio_obj_rep.each{|g_n,value|
 chr = value.primary_accession.split(":")[2]
 
 value.features.each {|featu|
  if featu.feature == "CTTCTT_match_chr"
    if featu.position.match(/complement/)
     coordenates = featu.position.split("(") #in case that the match is in the complement reverse strand
     coordenates = coordenates[1].chomp(")")
     start_posit,end_posit = coordenates.split("..")
    
    else 
     start_posit,end_posit = featu.position.split("..")
    end
    gff3_file << "Chr#{chr}\tBioRuby\trepeated_sequence\t#{start_posit}\t#{end_posit}\t.\t#{featu.assoc["strand"]}\t.\t#{featu.assoc["id"]}\n\n"
  end}
  
  }
gff3_file.close

