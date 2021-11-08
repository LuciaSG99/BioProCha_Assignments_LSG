## main code
file = ARGV[2]
require 'rest-client'
require '/home/osboxes/BioinformaticsCourseGitHub/BioProCha_Assignments_LSG/Assignment_2/InteractionNetwork_obj.rb'
require '/home/osboxes/BioinformaticsCourseGitHub/BioProCha_Assignments_LSG/Assignment_2/gene_obj.rb'



#saving the genes of the list in an array:
def read_file_genes(file)
  genes_list = Array.new
  lines = IO.readlines(file)
  lines.each {|gene| genes_list << gene.chomp}
  
return genes_list
end


#### OBTAIN THE INTERACTIONS OF THE GENES FROM THE ORIGINAL LIST
genes_list = read_file_genes(file)
hash_interactions = Hash.new
genes_list.each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 1)} #each value of the hash is a Gene object that correspond with the genes from the original list
 ### AQUI IMPRIME BIEN LOS GENES OBJETOS GEN CON TUS INTERACCIONES (en minus) Y SU NIVEL (int)
 #hash_interactions.values.each{|gene_original|
 # puts "Gene_object:"
 # puts gene_original.level
 # puts }

puts "Hash interaction length with only original genes:"
puts hash_interactions.length
#### STORE THE NEW GENES (Level 2) IN ANOTHER LIST, TRANSFORM THEM INTO A GENE OBJECT AND LOOK FOR THEIR INTERACTIONS
genes_list2 = Array.new()
hash_interactions.values.each {|gene_obj| gene_obj.interactions.each {|gene| genes_list2 << gene}}


genes_list2.each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 2)} #each value of the hash is a Gene object that correspond with the genes from the original list

#### TAMBIEN LO ALMACENA CORRECTAMENTE LOS GENES DEL SEGUNDO NIVEL CON SUS RESPECTIVAS INTERACCIONES Y NO CAMBIA EL NIVEL
#hash_interactions.values[50..89].each{|gene_original|
#  puts "Gene_object:"
#  puts gene_original.inspect
#    puts }
#puts hash_interactions.length
#puts "Hash interaction length with only second level genes:"
#puts hash_interactions.length

#### STORE THE NEW GENES (Level 3) IN ANOTHER LIST, TRANSFORM THEM INTO A GENE OBJECT AND LOOK FOR THEIR INTERACTIONS
genes_list3 = Array.new()
hash_interactions.values[168..-1].each {|gene_obj|
  if gene_obj.level == 2
    gene_obj.interactions.each {|gene| genes_list3 << gene}
  else
    next
  end}

genes_list3.each {|gene|
   if hash_interactions.has_key?(gene)
       next   
   else 
     hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 3)
   end} #each value of the hash is a Gene object that correspond with the genes from the original list

#hash_interactions.values[560..-1].each{|gene_original|
#  puts "Gene_object:"
#  puts gene_original.inspect
#  puts }
#puts "Length of hash:"
#puts hash_interactions.length

#puts "Hash interaction length with all genes:"
#puts hash_interactions.length
#### TAMBIEN LO ALMACENA CORRECTAMENTE LOS GENES DEL SEGUNDO NIVEL CON SUS RESPECTIVAS INTERACCIONES

### DELETING GENES THAT DON'T HAVE INTERACTIONS:
######## ESTO SUPUESTAMENTE TAMBIEN LO HACE BIEN, ELIMINA GENES SIN INTERACCIONES Y LOS GENES 2 Y 3 CON MENOS DE DOS INTERACCIONES 
hash_interactions_updated = Hash.new #new hash that will contain only those genes that have interactions
hash_interactions.values.each{|gene_obj|
  if gene_obj.interactions.length < 1
    next
  end
  if gene_obj.level > 1 and gene_obj.interactions.length < 2
    next
  end
  hash_interactions_updated[gene_obj.agi_locus] = gene_obj}

#hash_interactions_updated.values[380..-1].each{|gene_original| #to see if it worked
#   puts "Interactions :"
#  puts gene_original.inspect
#  puts }
##
#puts "Length of the hash updated: "
#puts hash_interactions_updated.length

######### CREATING AND OBTAINING INTERACTIONS NETWORKS USING THE INTERACTION NETWORK CLASS:


list_of_networks = Array.new
hash_interactions_updated.values.each {|gene|                            
                            if list_of_networks.any?{|int_network| int_network.network.any?{|int| gene.agi_locus == int.agi_locus}} #para no hacer redes repetidas
                              next
                            else
                              network = InteractionNetwork.new(:gene_original => gene,:all_genes => hash_interactions_updated)                               
                              list_of_networks << network
                            end}
#puts "FIRST NETWORK:"
#list_of_networks[0].network.each {|gene_obj|
# if gene_obj.level.to_i == 1
#   puts gene_obj.agi_locus
# end}
#puts
#
#puts "SECOND NETWORK:"
#list_of_networks[11].network.each {|gene_obj|
# if gene_obj.level.to_i == 1
#   puts gene_obj.agi_locus
# end}


#puts "Number of networks: "
#puts list_of_networks.length #--> cuantas redes obtengo
#list_of_networks.each{|red|
#     n = 0
#     puts "network: "
#     red.network.each{|int| if int.level.to_i == 1 
#       n+=1#ver cuantos genes originales tenemos en cada network conseguida 
#     end}
#     puts "number of original genes: " + n.to_s + "\t" + "Length of the network: " + red.network.length.to_s} ## ver la longitud de cada red

#### A FUNCTION THAT STORE IN AN ARRAY ONLY THE NETWORKS THAT INCLUDES AT LEAST TWO OF THE ORIGINAL GENES:
def remove_networks_no_functional(array_of_networks,new_list_networks)
  array_of_networks.each {|int_network_obj| int_network = int_network_obj.network 
    if int_network.count{|gene| gene.level.to_i == 1} < 2
       next
    else
       new_list_networks << int_network_obj
    end}
  return new_list_networks
end

list_networks_updated = remove_networks_no_functional(list_of_networks,Array.new) #applies the function to the list of total networks obtained
##
puts list_networks_updated.length #to see how many networks are obtained : ME APARECEN 2 REDES "LIMPIAS"
list_networks_updated.each{|element| puts "network:"
 puts element.kegg_annotations
 puts}}
  
  
