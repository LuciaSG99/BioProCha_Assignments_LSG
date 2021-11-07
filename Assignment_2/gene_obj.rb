class Gene
    attr_accessor :interactions
    attr_accessor :agi_locus
    attr_accessor :level
    require 'rest-client'
    
    def initialize(params = {})
        @agi_locus = params.fetch(:agi_locus,nil)
        @interactions = obtain_interactions(@agi_locus)
        @level = params.fetch(:level,nil).to_i #atribute to distinguish the origin of the gen. If it is from the original list = 1,
        #if it is a gene from the list_2 then = 2 and so on.  
    end
    
    def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
      response = RestClient::Request.execute({
        method: :get,
        url: url.to_s,
        user: user,
        password: pass,
        headers: headers})
      return response
      
      rescue RestClient::ExceptionWithResponse => e
        $stderr.puts e.inspect
        response = false
        return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
      rescue RestClient::Exception => e
        $stderr.puts e.inspect
        response = false
        return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
      rescue Exception => e
        $stderr.puts e.inspect
        response = false
        return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
    end 

    def obtain_interactions(agi_locus)
      result = fetch("www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{agi_locus}?format=tab25")  
      list = Array.new
      if result
          res = Array.new(result.body.split("\n"))
          res.each{|line| line = line.split("\t")
          if line[9].match(/Arabidopsis thaliana/) == nil
            next
          end
          if line[10].match(/Arabidopsis thaliana/) == nil
            next
          end          
          if line[4].match(/A[tT][1-5]g/) == nil
            next
          end
          if line[5].match(/A[tT][1-5]g/) == nil
            next
          end
          if line[14].match(/\d.\d{2}/)[0].to_f < 0.45
            next
          end    
          int1 = line[4].match(/A[tT][1-5]g\d+/)[0]
          
          int2 = line[5].match(/A[tT][1-5]g\d+/)[0]
          
          if list.include? int1 #to avoid insert repeated interactions
            next
          end         
          
          if list.include? int2
            next
          end
          
          if int1.casecmp(@agi_locus) == 0 #to ignore the case of the letters. If the strings are equals, it returns 0
             list << int2.downcase
          else
              list << int1.downcase
          end}
          
      else    
        list << ""
      end
         return list.select{|int| int != agi_locus}       
    end
    
end