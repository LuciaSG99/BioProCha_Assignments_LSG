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

genes_list = read_file_genes(file)
hash_interactions = Hash.new(0)

#### OBTAIN THE INTERACTIONS OF THE GENES FROM THE ORIGINAL LIST
genes_list[0..30].each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 1)} #each value of the hash is a Gene object that correspond with the genes from the original list
#lo pongo en minusculas las keys para quitarme los problemas de mezcla mayusculas y minusculas

#### STORE THE NEW GENES (Level 2) IN ANOTHER LIST, TRANSFORM THEM INTO A GENE OBJECT AND LOOK FOR THEIR INTERACTIONS
genes_list2 = Array.new()
hash_interactions.values.each {|gene_obj| gene_obj.interactions.each {|gene| genes_list2 << gene}}

genes_list2.each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 2)} #each value of the hash is a Gene object that correspond with the genes from the original list

#### STORE THE NEW GENES (Level 3) IN ANOTHER LIST, TRANSFORM THEM INTO A GENE OBJECT AND LOOK FOR THEIR INTERACTIONS
genes_list3 = Array.new()
hash_interactions.values.each {|gene_obj|
  if gene_obj.level == 2
    gene_obj.interactions.each {|gene| genes_list3 << gene}
  else
    next
  end}

genes_list3.each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 3)} #each value of the hash is a Gene object that correspond with the genes from the original list


######### CREATING AND OBTAINING INTERACTIONS NETWORKS USING THE INTERACTION NETWORK CLASS:
list_of_networks = Array.new
hash_interactions.values.each {|gene| network = InteractionNetwork.new(:gene_original => gene,:all_genes => hash_interactions)
                                    list_of_networks << network}

# puts list_of_networks.length #--> ME SALEN 21 REDES SIN LIMPIAR
 list_of_networks.each{|red|puts "network:"
   puts red.network   
   puts "\n"} ## PARA VER COMO SON ESAS REDES SIN LIMPIAR

#### A FUNCTION THAT STORE IN AN ARRAY ONLY THE NETWORKS THAT INCLUDES AT LEAST TWO OF THE ORIGINAL GENES:
def remove_networks_no_functional(array_of_networks,list_original_genes,new_list_networks)
  array_of_networks.each {|int_network_obj| int_network = int_network_obj.network 
    if int_network.count{|gene| list_original_genes.include?(gene)} < 2
       next
    else
       new_list_networks << int_network_obj
    end}
  return new_list_networks
end

#original_genes = genes_list.each {|gene| gene.downcase}
#list_of_int_networks = remove_networks_no_functional(list_of_networks,original_genes,Array.new) #applies the function to the list of total networks obtained
##
#puts list_of_int_networks.length #to see how many networks are obtained : ME APARECEN 0 REDES "LIMPIAS"
#list_of_int_networks.each{|element| element.network #to see what contain each network
#  puts}
#  
  
# hash = { "AT1G23094": Gene_obj,...}