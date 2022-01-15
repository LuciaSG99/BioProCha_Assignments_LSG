### ASSIGNMENT 3 : MAIN CODE
require 'rest-client'  
require 'bio'
list_genes = ARGV[0]

seq_positive = Bio::Sequence.auto("GAGAAGAATAGCTGTCTTCTTGG")



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
genes[0..4].each {|gene| array = Array.new
  array << create_obj_embl(gene)
  hash_embl_objs[gene.upcase] = array}


########## 2. SEARCH FOR THE CTTCTT SEQUENCE IN EACH EXON:

#Store the exons sequences for each gene:
hash_embl_objs.each{|key, list|
  re = Regexp.new(key)
  hash_exon = Hash.new #hash that will contain all the exons of each gene with its respective sequence 
  list[0].ft.each {|feature|
    if feature.feature == "exon" 
      if feature.qualifiers[0].value.match(re)
        if feature.position.match(/complement\(\d+\.\.\d+\)/) == nil  and feature.position.match(/\d+\.\.\d+/) == nil #ŧo avoid repeated genes or some other elements 
            next
        end
        hash_exon[feature.qualifiers[0].value.match(/exon\d/)[0]] = list[0].seq.splice(feature.position) #splice give us the part of the original sequence where the exon is
      end
    end}
  hash_embl_objs[key] << hash_exon} 

#Seach for the repetitions in the exons:
regular_expression = Bio::Sequence::NA.new("CTTCTT") 

def find_seq_exon(regexp,seq)
  re = Regexp.new(regexp.to_re)
  match = seq.match(re)
  if match == nil
    return   
  else
  return match.offset(0)
  end
end

#puts hash_embl_objs.values
hash_embl_objs.each {|key,list| #puts "gene:"+ key.to_s
 # puts "exons:"
  hash_seq_exons = list[1]
  hash_seq_exons.each {|key,exon| position_match = find_seq_exon(regular_expression,exon) 
  #puts key
  #puts position_match
  }}
puts "coordenates of the exon inside the sequence: [2..17]"
puts "position match inside all sequence:"
puts seq_reverse
puts find_seq_exon(regular_expression,seq_reverse)
puts "position match inside the exon:"
prueba1 = seq_reverse.seq.splice("2..17")
puts prueba1
prueba = find_seq_exon(regular_expression,prueba1) 
puts prueba
#puts hash_embl_objs.values[1][1] #lista con las posiciones de los exones (de todos!!) por cada gen
#puts hash_embl_objs.keys[1]

#puts hash_embl_objs.values[1][1]
#lista con la posición del exon 
#contemplar la idea de en vez de meter las posiciones como una lista, hacer un diccionario donde cada llave es un exon que contiene una lista donde el primer elemento
#es la posicion de start del exon y el segundo elemento es "end"

feature_exon = Bio::Feature.new(feature = "CTTCTT_match",position = coordenates_forward)
              feature_exon.append(Bio::Feature::Qualifier.new('strand', '+')) #indicates that it is in the forward strand
              feature_exon.append(Bio::Feature::Qualifier.new('id', "#{key}_#{feature.qualifiers[0].value.match(/exon\d/)[0]}")) #GeneName_Exon as an id
              bioseq_obj = embl_obj.to_biosequence #transforms the EMBL object into a Bio sequence object
              bioseq_obj << feature_exon}
              
              
if feature.feature == "CTTCTT_match"
    #puts feature.inspect
     puts "position of CTTCTT in #{feature.assoc["id"]} is #{feature.position}"
     puts "\n"
  end}



def obtain_position_cttctt(gene,embl_obj,exon,chr_start)
  re = Regexp.new(key)
  position_matches = Array.new() #Array that will contain the coordenates of all the matches for one gene
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
        seq_exon = embl_obj.seq.splice(position) #store the sequence of the exon
        
        if localization_obj.strand == 1 #if the exon is in the forward strand            
                      
            matches = seq_exon.seq.to_enum(:scan,/(?=(cttctt))/).map{Regexp.last_match} #to get all the matches in overlapping situations like this: CTTCTTCTT 
            
            unless matches.empty? #unless no matches were found:      
               
              matches.each{|match|
              
            ###Obtain the coordenates of the match in the whole gene sequence:
            #to obtain the coordenates of the match in the whole sequence, we have to sum the initial position (nucleotide) of the exon plus the initial
            #position of the match.
              if exon == 1
                start_pos_seq = match.begin(0).to_i + localization_obj.from.to_i #initial position of the match in the sequence (GENE)
              
              else
                start_pos_seq = match.begin(0).to_i + chr_start.to_i #initial position of the match in the whole chromosome sequence 
              end
                                           
              end_pos_seq = start_pos_seq.to_i + 5 #the final position of the match is calculated by adding 5 to the start position (because the repeat has 6 nucleotides)
              coordenates = "#{start_pos_seq.to_s}..#{end_pos_seq.to_s}"
              position_matches << coordenates
              
            end
        
        end  
            
        if localization_obj.strand == -1 #if the exon is in the complementary reverse strand
                    
            #Because the exon is in the complementary reverse strand, we can obtain the position of the match CTTCTT by searching for the sequence "AAGAAG" in
            #the forward strand. 
           matches_2 = seq_exon.seq.to_enum(:scan,/(?=(aagaag))/).map{Regexp.last_match}  
                      
           unless matches_2.empty? #unless no matches were found:      
                      
             matches_2.each{|match|
             if exon == 1
                start_pos_seq = match.begin(0).to_i + localization_obj.from.to_i #initial position of the match in the sequence (GENE)
              
              else
                start_pos_seq = match.begin(0).to_i + chr_start.to_i #initial position of the match in the whole chromosome sequence 
              end
                                           
              end_pos_seq = start_pos_seq.to_i + 5 #the final position of the match is calculated by adding 5 to the start position (because the repeat has 6 nucleotides)
              coordenate = "complement(#{start_pos_seq.to_s}..#{end_pos_seq.to_s})"
              position_matches << coordenates
              
  return position_matches            
end