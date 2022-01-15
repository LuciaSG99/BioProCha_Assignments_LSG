list_networks_updated.each{|element| puts "network:"
 puts element.kegg_annotations
 puts
 }
    
###### CREATING REPORT:
report = File.open(seed_stock_modified,"w")
report << "REPORT: IDENTIFICATION OF INTERACTION NETWORKS\n\n" #header
report << "#{list_networks_updated.length} networks have been found.\n"
list_networks_updated.each {|net| report << "\n\nNetwork: Contains #{net.network.length} genes in total\n"
report << "The original genes that participe in this network are: \n"
net.each{|gene| if gene.level.to_i == 1
  report << gene.agi_locus.upcase-to_s + "\n"}
report << "\n" + "The KEGG annotations of the original genes involved in this network are:\n"
net.kegg_annotations.each {|key,value| report << "KEGG ID: " + key.to_s + "\t" + "Pathway name: " + value.to_s}
report << "\n" + "The GO annotations of the original genes involved in this network are:\n"
net.go_annotations.each {|key,value| report << "GO ID: " + key.to_s + "\t" + "Pathway name: " + value.to_s}}

#puts "KEGG ANNOTATIONS: "
#list_networks_updated.each{|net| puts net.kegg_annotations
# puts
# puts
# }
#puts "GO ANNOTATIONS: "
#list_networks_updated.each{|net| puts net.go_annotations
# puts
# puts
# } 
 #puts "Number of networks: "
#puts list_networks_updated.length #--> cuantas redes obtengo
#list_networks_updated.each{|red|
#     n = 0
#     puts "network: "
#     red.network.each{|int| if int.level.to_i == 1 
#       puts int.agi_locus
#       n+=1#ver cuantos genes originales tenemos en cada network conseguida 
#     end}
#     puts "number of original genes: " + n.to_s + "\t" + "Length of the network: " + red.network.length.to_s
#     puts "\n"} ## ver la longitud de cada red
    hash_updated_sorted.values.each {|gene|                            
                            if list_of_networks.each{|int_network| int_network.network.any?{|int| int.agi_locus.to_s == gene.agi_locus.to_s}} #para no hacer redes repetidas
                              next
                            else
                              network = InteractionNetwork.new(:gene_original => gene,:all_genes => hash_updated_sorted)                               
                              list_of_networks << network
                            end}
    
    list_networks_updated = remove_networks_no_functional(list_of_networks,Array.new) #applies the function to the list of total networks obtained
puts "Number of networks: "
puts list_networks_updated.length #--> cuantas redes obtengo
list_networks_updated.each{|red|
     n = 0
     puts "network: "
     red.network.each{|int| if int.level.to_i == 1 
       puts int.agi_locus
       n+=1#ver cuantos genes originales tenemos en cada network conseguida 
     end}
     puts "number of original genes: " + n.to_s + "\t" + "Length of the network: " + red.network.length.to_s
     puts "\n"} ## ver la longitud de cada red
hash_updated_sorted.values.each {|gene| network = InteractionNetwork.new(:gene_original => gene,:all_genes => hash_updated_sorted)                           
                            if network.network.any?{|int| int.agi_locus.to_s == gene.agi_locus.to_s} #para no hacer redes repetidas
                              next
                            else                                                        
                              list_of_networks << network
                            end}