### ASSIGNMENT 3 : MAIN CODE
require 'rest-client'  
require 'bio'
list_genes = ARGV[0]

######## 1. RETRIEVE THE SEQUENCES OF THE GENES: 

#Function that saves the genes of the list in an array:
def read_file_genes(file)
  genes_list = Array.new
  lines = IO.readlines(file)
  lines.each {|gene| genes_list << gene.chomp}
  return genes_list
end

genes = read_file_genes(list_genes)

#Function that creates a EMBL object 
def create_obj_embl(gene)
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}"  # create a "URI" object (Uniform Resource Identifier: https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)
  response = RestClient.get(address)
  record = response.body
  entry = Bio::EMBL.new(record)
  return entry
end

#Creation of a Hash which will contain the genes of the list as a Embl objects from BioRuby, in the method ".seq" we can retrieve the sequence:
hash_embl_objs = Hash.new
genes[0..8].each {|gene| hash_embl_objs[gene.upcase] = create_obj_embl(gene)}


########## 2 AND 3. SEARCH FOR THE CTTCTT SEQUENCE IN EACH EXON AND CREATE NEW SEQUENCE FEATURE:

#FUNCTION TO CREATE NEW FEATURES : 
def new_feature(name,coordenates,strand,id_gene_exon)
  feature = Bio::Feature.new(feature = name,position = coordenates)
  feature.append(Bio::Feature::Qualifier.new('strand', strand)) #indicates that it is in the forward strand
  feature.append(Bio::Feature::Qualifier.new('id', id_gene_exon)) #GeneName_Exon as an id
  return feature             
end

########Store the exons sequences for each gene and seach for the repetitions in the exons:
hash_bio_seq_obj = Hash.new #to add the new features we have to use Bio::Sequence objects instead of Bio:EMBL objects. This hash will contain those objects 
hash_embl_objs.each{|key, embl_obj| #for each key (Gene name, for example = "AT4g27030") and list containing 
  
  bioseq_obj = embl_obj.to_biosequence #transforms the EMBL object into a Bio sequence object
  
  re = Regexp.new(key)
  embl_obj.ft.each {|feature|
    if feature.feature == "exon" #if the feature is a exon:
     
      if feature.qualifiers[0].value.match(re)#if it is a exon that belongs to the gene of the list:
      # puts feature.qualifiers[0].value.match(re)[0]
        if feature.position.match(/complement\(\d+\.\.\d+\)/) == nil  and feature.position.match(/\d+\.\.\d+/) == nil #ŧo avoid repeated genes or some other elements 
            next
        end
        localization_obj = feature.locations.locations[0]
        start = localization_obj.from.to_s
        pos_end = localization_obj.to.to_s
        
        position = start + ".." + pos_end #coordenates of the exon
        
        if localization_obj.strand == 1 #if the exon is in the forward strand
            seq_exon = embl_obj.seq.splice(position) #store the sequence of the exon 
                      
            matches = seq_exon.seq.to_enum(:scan,/(?=(cttctt))/).map{Regexp.last_match} #to get all the matches in overlapping situations like this: CTTCTTCTT 
            
            unless matches.empty? #unless no matches were found:      
               
              matches.each{|match|
              
            ###Obtain the coordenates of the match in the whole gene sequence:
            #to obtain the coordenates of the match in the whole sequence, we have to sum the initial position (nucleotide) of the exon plus the initial
            #position of the match.
              
              start_pos_seq = match.begin(0).to_i + localization_obj.from.to_i #initial position of the match in the sequence 
               
              end_pos_seq = start_pos_seq.to_i + 5 #the final position of the match is calculated by adding 5 to the start position (because the repeat has 6 nucleotides)
              coordenates_forward = "#{start_pos_seq.to_s}..#{end_pos_seq.to_s}"
              
              #creates the features object
              feature_match = new_feature("CTTCTT_match",coordenates_forward,'+',"#{key}_#{feature.qualifiers[0].value.match(/exon\d/)[0]}")
              
              #if ft1.position == feature_match.position && feature_match.qualifiers[1].qualifier == ft1.qualifiers[1].qualifier
               
               bioseq_obj.features << feature_match
              #end
             }
            end
        
        end  
            
        if localization_obj.strand == -1 #if the exon is in the complementary reverse strand
           seq_exon = embl_obj.seq.splice(position) #store the sequence of the exon (in the forward strand!)
           
            #Because the exon is in the complementary reverse strand, we can obtain the position of the match CTTCTT by searching for the sequence "AAGAAG" in
            #the forward strand. 
           matches_2 = seq_exon.seq.to_enum(:scan,/(?=(aagaag))/).map{Regexp.last_match}  
                      
           unless matches_2.empty? #unless no matches were found:      
                      
             matches_2.each{|match|
                            
           ###Obtain the coordenates of the match in the whole gene sequence:
             
             start_pos_seq = match.begin(0).to_i + localization_obj.from.to_i #initial position of the match in the sequence --> LO HACE BIEN
             
             end_pos_seq = start_pos_seq.to_i + 5 #the final position of the match is calculated by adding 5 to the start position (because the repeat has 6 nucleotides)
             coordenates_forward = "complement(#{start_pos_seq.to_s}..#{end_pos_seq.to_s})" #ESTO TB
            
             #creates the features object
             feature_match = new_feature("CTTCTT_match",coordenates_forward,'-',"#{key}_#{feature.qualifiers[0].value.match(/exon\d/)[0]}")
             
             #to avoid repeating features:
             
             
             # if ft1.position == feature_match.position && feature_match.qualifiers[1].qualifier == ft1.qualifiers[1].qualifier
               
               bioseq_obj.features << feature_match
             # end
             }
           end
        end
             
      end
    end}
  hash_bio_seq_obj[key] = bioseq_obj} 


###### Search for genes that do not contain the sequence in their exons:
list_genes_no_repeat = Array.new #Array that will contain those genes

 hash_bio_seq_obj.each{|key,bio_obj| #for each gene and its Biosequence object:
  #if that gene does not have features called CTTCTT_match, it means that do not have that sequence
  #in its exons because that type of feature was only created when there was a match
  if bio_obj.features.any?{|feature| feature.feature == "CTTCTT_match"} == false 
   list_genes_no_repeat << key #store the name of the gene in the array
  end}

#### Store only the genes that have the sequence CTTCTT at least in one: 
hash_bio_obj_rep = Hash.new
hash_bio_seq_obj.each{|key,obj|
 unless list_genes_no_repeat.include? key
  hash_bio_obj_rep[key] = obj 
 end}

 #puts hash_bio_obj_rep.values[0].features.inspect
 
 
########### 4A. CREATING A GFF3-FORMATTED FILE OF THE FEATURES:

gff3_file = File.open("new_features.gff", 'w')
gff3_file << "##gff-version 3\n\n"
hash_bio_obj_rep.each{|key,value|
 value.features.each {|featu|
  if featu.feature == "CTTCTT_match"
    if featu.position.match(/complement/)
     coordenates = featu.position.chomp("complement(") #in case that the match is in the complement reverse strand
     coordenates = coordenates.chomp(")")
     start_posit,end_posit = coordenates.split("..")
    
    else 
     start_posit,end_posit = featu.position.split("..")
    end
    gff3_file << "#{key}\tBioRuby\trepeated_sequence\t#{start_posit}\t#{end_posit}\t.\t#{featu.assoc["strand"]}\t.\t#{featu.assoc["id"]}\n\n"
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

hash_embl_objs




