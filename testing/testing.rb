
#Store the exons coordenates for each gene:
hash_embl_objs.each{|key, list|
  re = Regexp.new(key)
  #array_positions_exon = Array.new

  hash_position_exon = Hash.new
  list[0].ft.each {|feature|
    if feature.feature == "exon"
      
      if feature.qualifiers[0].value.match(re)
        #puts feature.qualifiers[0].value.match(re)[0]
        positions_exon = feature.position.split("..") #array that contains the start and end position of the end as the first and second element respect
        hash_position_exon[feature.qualifiers[0].value.match(/exon\d/)[0]] = positions_exon
        #array_positions_exon << feature.position
      end
    end}
  hash_embl_objs[key] << hash_position_exon} 

def obtain_seq_exons(seq,exon_list)
  
    
end

#Seach for the repetitions in the exons:
regular_expression = Bio::Sequence::NA.new("CTTCTT") 

def find_seq_exon(regexp,seq)
  re = Regexp.new(regexp.to_re)
  match = seq.match(re)
  return match
end

hash_embl_objs.values.each {|list| obj_embl = list[0]
seq_exons = obj_embl.seq
matches_seq = find_seq_exon(regular_expression,seq_exons)


#hash_embl_objs.values.each{|obj_embl| puts "Sequence: #{obj_embl.entry_id}"
  #puts "+ strand:"
  #puts obj_embl.seq
  #puts "\n - strand:"
  #puts obj_embl.seq.complement
 # puts "\n\n\n"}