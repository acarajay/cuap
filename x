>JX275887.1 Cocos nucifera beta-ketoacyl-ACP synthase I (KASI) mRNA, complete cds
ATCCCCGTATCTTCTCACCTCCTCCTCTCTTGCTTATATCCTGCACTCTCCCCCTCTCCACACCTCTCCA
>DQ861409.1 Cocos nucifera WRKY2 gene, partial cds
GGTTGTGCTCGCCCTCGTACGTTACGATCAGCATCGACGGATCGTCGGGCGCCCGCTCCACGTGCTTCC
>DQ861408.1 Cocos nucifera WRKY1 gene, partial cds
GTTGTGCTCGCCCTC

JX275887.1, DQ861409.1, DQ861408.1 

>KU711779.1 Cocos nucifera oleosin isoform 500c mRNA, complete cds
ATGGCGGAGCAGCAGCCAACGTCGCACAGGGTGGTGAAGGGGGTGACGGCGGCGACCATTGGCGGATCGC
TGCTGCTGCTCTCGGGGCTAACGCTGGCGGGGACGGTGATCGGGCTGGCGGTGGTGACACCGCTGCTGGT
CATCTTCAGCCCGGTGCTGGTGCCGGCGGTCATCACCGTGTTCCTCCTGGTGACGGGCTTCGTCACCTCC
GGCGTGTTAGGGGTGGCGGCGCTCTCGGTGCTGTCCTGGTTGTACAAGTACCTCACCGGCAAGCGCGTAC
CGGGGGCCGAGCAGCTGGAACAAGCACGGGCCCGCTTCGCCTCCAAAGCCCGCGACGTCAAGGAGTCCGC
TCAGCACCGCACCGAGCAGGCCTAG

>KU711778.1 Cocos nucifera oleosin isoform 500a mRNA, complete cds
ATGGCGGAGCAGCAGAAGGAGCGCCTGGCGACGTCGCACGCAGTTGTGAGGGGAGTGATGGCGGCGACCA
TCGGGGGATCGCTGCTGCTGCTGTCGGGGCTAACGCTGACTGGCACGGTAATCGGGCTGACGGTGCTGAC
GCCGCTGCTGGTCCTCTTCAGCCCGGTGCTGGTGCCGGCGGCCATCGCCGTGTTCCTGCTGGTGGCCGGC
TTCGTCACCTCCGGCGGGTTCGGGTTGGCGGCGCTATCGGTGCTCTCCTGGATGTACAAGTACCTCACCG
GCAGGCGCCCGCCGGGGTCCGAGCAGCTGGAGCAGGCACGCGCCCGGCTTGCCTCCAAGGCCCGCGACAT
CAAGGAGTCCGCGCAGCACCGCATCGACCAGGCCCAGTCCTCCTAA




















 
 
 	    #---------------Set the input values as reactive ------------------
	    datasetInput <- reactive({
	    	one <- input$inputfile
	    	#two <- input$two
	    	#three <- input$three


  		})

  		output$summary <- renderPrint({
		    dataset <- datasetInput()
			if(is.null(dataset)){return()}
			print(is.null(dataset))
			print(dataset)
		    #summary(dataset)
		})

			    #This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file
	    #file$datapath -> gives the path of the file
	    data <- reactive({
	      file1 <- input$file
	      seq_coco <- read.GenBank(input$accessno)

	      if(is.null(file1)){return()}
	      read.table(file=file1$datapath, sep = "\t", header = FALSE, stringsAsFactors =FALSE)
	      
	    })


	    output$sample <- renderPrint({
		    genome <- input$sequence 

            if(is.null(genome)){return()}
	        #write input sequence into genome.fasta
	        writeGenome <- write.fasta(sequences = genome, names = names(genome), nbchar = 60, file.out = "genome.fasta")
	        
	        #read the input sequence from a genome.fasta, sequence as DNA, string, uppercase, sequence attributes should be set, only sequences as returned , the '>' at the beginning of the description lines is removed in the annotations of the sequences
	        dnafile <- read.fasta(file = "genome.fasta", seqtype = c("DNA"), as.string = TRUE, forceDNAtolower = FALSE, set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = TRUE)
	        #print(dnafile)
	        #print(length(dnafile))
	        #read each sequences in dnafile
	        #print(dnafile[2])
	        #print(getName(dnafile))
	        #print(getAnnot(dnafile[[2]]))
	        
			#create a new list containing only the sequences
	        seqlist=''
	        allcodons=''
	        allaminoacids=''
	        annotList=''
	        i = 2
	        j = 1
	        for(i in 2:length(dnafile)){
	          seqlist[j] <- c(dnafile[i])
	          #print(getAnnot(dnafile))
	          annotList[j] <- c(getAnnot(dnafile[[i]]))
	          #print(seqlist[j])
	          #print(annotList[j])
	          j = j+1
	        }
	        
	        #for each sequence in seqlist extract each codon
	        a = 1
	        #print(length(seqlist))
	        dnaletters=''
	        for(a in 1:length(seqlist)){   #index of each sequence in seqlist
	          #print(seqlist[a])  # but sequences here is of string value
	          #extract each character in a sequence of seqlist
	          invalid=FALSE
	          

	          seqstring <- trimSpace(seqlist[a])   # trims (" Some text. " == "Some text.")
	          newseq <- str_replace_all(seqstring, "[ \r\n]" , "")   #replace all newline, tabs and whitespace with ""
	          #print(newseq)
	          bases <- unlist(str_split(newseq, ""))  #[[1]] splits into individual letters/characters
	          dnaletters <- toupper(bases)
	          #print(dnaletters)
	          #newy <- gsub("[[:space:]]", "", dnaletters)
	          #if()
	          #print(newy)
	          
	          #---------prints annotation of the sequence-----------
	          print(annotList[[a]])
	          #-------------------------------------------
	          #extract each codon for each sequence in seqlist
	          #print(length(dnaletters))
	          b = 1
	          codon=''
	          aminoacid=''
	          codonlist <- hash("ATA"="I","ATC"="I","ATT"="I","ATG"="M",
	                            "ACA"="T","ACC"="T","ACG"="T","ACT"="T",
	                            "AAC"="N","AAT"="N","AAA"="K","AAG"="K",
	                            "AGC"="S","AGT"="S","AGA"="R","AGG"="R",
	                            "CTA"="L","CTC"="L","CTG"="L","CTT"="L",
	                            "CCA"="P","CCC"="P","CCG"="P","CCT"="P",
	                            "CAC"="H","CAT"="H","CAA"="Q","CAG"="Q",
	                            "CGA"="R","CGC"="R","CGG"="R","CGT"="R",
	                            "GTA"="V","GTC"="V","GTG"="V","GTT"="V",
	                            "GCA"="A","GCC"="A","GCG"="A","GCT"="A",
	                            "GAC"="D","GAT"="D","GAA"="E","GAG"="E",
	                            "GGA"="G","GGC"="G","GGG"="G","GGT"="G",
	                            "TCA"="S","TCC"="S","TCG"="S","TCT"="S",
	                            "TTC"="F","TTT"="F","TTA"="L","TTG"="L",
	                            "TAC"="Y","TAT"="Y","TAA"="_","TAG"="_",
	                            "TGC"="C","TGT"="C","TGA"="_","TGG"="W")
	
	          
	          for(b in 1:length(dnaletters)){
	          	#print(dnaletters[b])    # "A","T","G","C"
	          	if(dnaletters[b]=="A" || dnaletters[b]=="T" || dnaletters[b]=="G" || dnaletters[b]=="C"){
	          		if(b%%3==0){   # every 3rd character
		              nucleotide <- c(dnaletters[(b-2):b])    #"A""T""G"
		              y <- paste(nucleotide, collapse="")  #"ATG"
		              codon <- append(codon, y)
		              aminoacid <- append(aminoacid, codonlist[[y]])              #prints list of protein ("M","M""V","M""V""S")
		              
		            }	
	          	}else{
	          		print(c("Sequence contains invalid input:", dnaletters[b]))
	          		invalid = TRUE
	          		break
	          	}	            
	          }
	          

	          if(invalid==TRUE){
	          	break
	          }else{        
		        #solve for the codon count 
		          
		        proteinseq <- paste(aminoacid, collapse="")
		        if(codon[1]=="" || protein[1]==""){
		        	codon <- codon[2:length(codon)]
		        	aminoacid <- aminoacid[2:length(aminoacid)]
		        }
		        #print(listofcodons)  
		        #print(nucleotide)
		        #print(length(nucleotide))
		        #print(proteinseq)
		        #print(length(aminoacid))
		        tb <- data.frame(codon,aminoacid)

		        #View(tb)
		        #print(tb)
		        codoncount <- count(tb, c('codon','aminoacid'))   #counts for the frequency
		        print(codoncount)

		        #----------write into csv file the frequency counts ------------
		        write.table(codoncount, file="widefinaldata.csv", sep=",", row.names=FALSE)
		        ds1<-read.csv("widefinaldata.csv",header=TRUE, sep=",", row.names=1)

		        #ala<-ds2[,1:4]
		        #print(ala)	          	
	          }
	          
	          allcodons <- append(allcodons,codon)  #appends each codons per sequence in a given input
	          allaminoacids <- append(allaminoacids,aminoacid)
	        }
	        CODONS <- allcodons[2:length(allcodons)]
	        AMINO_ACID <- allaminoacids[2:length(allaminoacids)]
	        #print(CODONS)   #prints all codons in a given input except null values
	        cut <- data.frame(CODONS,AMINO_ACID)
	        newcut <- count(cut,c('CODONS', 'AMINO_ACID'))
	        newcut <- count(cut,c('CODONS','AMINO_ACID'))
	        #View(newcut)
	        wds<- as.data.frame(spread(newcut, CODONS, freq, fill=0))     #
			View(wds)
	        write.table(wds, file="codon_usage_table.csv", sep=",", row.names=FALSE)
			ds2<-read.csv("codon_usage_table.csv",header=TRUE, sep=",", row.names=1)


            ds3 <- data.matrix(ds2)
            View(ds3)
			#hm <- heatmap.2(ds3, main="Sample", trace="none", margins=c(10,12))
			#print(hm)

			#hm = heatmap(ds3,scale="column",col=heat.colors(256),main="Codon Usage Table", Rowv=NA,Colv=NA)
			#hm
			ala <- ds3[,13:15]
			cys <- ds3[2,]
			asp <- ds3[3,]
			#function for generating heatmaps for the codon usage
			heatmap<-function(data) {
			    my_data<-data
				res.dist<-get_dist(my_data,stand=TRUE,method="spearman")
				#Visualize the dynamics of the distances
				fviz_dist(res.dist,gradient=list(low="#00AFBB",mid="white",high="#FC4E07"))
				#High "Orange"; LOw "Blue
			}
			hm = heatmap(ds3)
			hm
			#heatmap(cys)
		

    	})
