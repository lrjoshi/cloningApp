library(shiny)
library(dplyr)
library(stringr)
library(DT)

relist<- read.csv("re.csv",header=T)



#first function for complement 
mycomplement <-function (d){
  #convert to upper characters
  d=toupper(d)
  #start empty variable for complement
  complement =as.character("")
  for (i in 1 :nchar(d)) {
    if (substring(d,i,i)=="A") {
      complement1 =as.character ("T")
      
    }else if (substring (d,i,i)=="G"){
      complement1 =as.character("C")
      
    }else if (substring (d,i,i)=="T"){
      complement1 =as.character("A")
      
    }else if (substring (d,i,i)=="C"){
      complement1 =as.character ("G")
    }
    complement=paste(complement,complement1)
  }
  
  remspace=gsub(" ","",complement)
  print(remspace,quote=FALSE)
  
}



#now make reverse  of complement sequence 
#first lets run the function for reverse complement 

reverse=function (d){
  complement=tolower(d)
  complement=gsub(" ","",complement)
  start1=as.character("")
  for (i in 0: nchar(complement)) {
    getletters=as.character(substring(complement,(nchar(complement)-i),(nchar(complement)-i)))
    start1=paste(start1,getletters)
  }
  #gsub command to remove spaces between the letters 
  print(gsub(" ", "", start1), quote=F)
}




  
ui <- fluidPage(theme = shinytheme("cyborg"),
  titlePanel("Shiny App for Designing Primers for DNA Cloning"),
  sidebarLayout(
    sidebarPanel (
      textInput("inputDNA",strong("Paste your DNA sequence below",style="color:green"),value="AAGGCCTTGGACCGGTTAGGACGCTGAAGCTTGGCCGGAATTTGGATGGACCACAGGATTAGACCAGGGAGGATCCGGAAGGATGAGGACCAGGATTAGAAGCTTGGAA"),
      br(),
      uiOutput("enzyme1Output"),
      uiOutput("enzyme2Output"),
      #selectInput("enzyme1","Enzyme for forward primer",c("HindIII-HF","KpnI-HF")),
      #selectInput("enzyme2","Enzyme for reverse primer",c("HindIII-HF","KpnI-HF")),
      numericInput("forward",strong("No of basepairs in forward primer",style="color:green"),22 , min=18, max=100),
      numericInput("reverse",strong("No of basepairs in reverse primer",style="color:green"),22, min=18,max=100)
      
      ),
    

   
      mainPanel(
       
        tableOutput("table"),
        hr(),
        br(),
        textOutput("site1"),
        tableOutput("search1"),
        
        
        br(),
        br(),
        textOutput("site2"),
        tableOutput("search2"),
        
       
        br(),
        hr(),
        strong("Cleavage sites of the selected enzymes",style = "color:blue"),
        tableOutput("cut_list"),
        
        hr(),
        h6("Powered by R Studio"),
        img(src="shiny.png",width=100,height=100,alt="Face"),
        h6("Developed by: Lok Raj Joshi",style="color:blue")
        )
      )
)
  



server <- function(input,output,session) {
  
  
  output$enzyme1Output <- renderUI({
    selectInput("enzyme1", strong("Select Enzyme for the Forward Primer",style="color:green"),
                sort(unique(relist$name)),
                selected = "HindIII-HF")
  })
  
  output$enzyme2Output <- renderUI({
    selectInput("enzyme2", strong("Select Enzyme for the Reverse Primer",style="color:green"),
                sort(unique(relist$name)),
                selected = "BamHI-HF")
  })
  
  

  
  
  
  
  
    
  
  
  
  
  
  
  output$table <- renderTable({  
  #required inputs
    enzyme1 = input$enzyme1
    enzyme2= input$enzyme2
    inputDNA= input$inputDNA
    inputDNA=toupper(inputDNA)
    forward= input$forward
    reverse=input$reverse
    length= nchar(inputDNA)
    library(dplyr)
    library(dplyr)
    
    #extract row containing enzymes
    e1=filter(relist, name == enzyme1)
    e2=filter(relist, name== enzyme2)
    
   
    #now extract restriction sites sequence 
    e1_sequence = as.character(tolower(e1$seq))
    e2_sequence = as.character (tolower(e2$seq))
    
    #extract forward primer sequence 
    library (stringr)
    fw= str_sub(inputDNA, start=1, end= forward)
    
    
    #extracting sequence from last 
    rv= str_sub (inputDNA, -reverse)
    
    
   
    
    #get complement sequence of the primer 
    complemement_sequence = mycomplement(rv)
    
    
  
    #now get the reverse complement 
    reverse_complement_sequence <-reverse(complemement_sequence)
    
    
    #finalizing forward primer , add 5 extra base pairs at the beginning 
    fw_primer= paste("AGGAC",e1_sequence,fw)
    
    #nice formatting
    fw_primer1=str_c (gsub(" ","",fw_primer), sep = "",collapse = NULL)
    fw_primer1=print(fw_primer1,quote=F)
    
  
    
    
    #get the reverve primer sequence 
    rev_primer= paste("ACAGC",e2_sequence,toupper(reverse_complement_sequence))
    
    
    #nice formatting
    rev_primer1=str_c (gsub(" ","",rev_primer), sep = "",collapse = NULL)
    rev_primer1=print(rev_primer1,quote=F)
    
    
 #Generate table    
  Forward=fw_primer1
  Reverse=rev_primer1
  table=data.frame(Forward,Reverse)
  table
  
})   
 
  
  
output$site1 <-renderText({
  paste("Location of", input$enzyme1, "enzyme within input DNA sequence")}) 
output$search1 <- renderTable({
  
  #required inputs
  enzyme1 = input$enzyme1
  enzyme2= input$enzyme2
  inputDNA= input$inputDNA
  inputDNA=toupper(inputDNA)
  forward= input$forward
  reverse=input$reverse
  length= nchar(inputDNA)
  library(dplyr)
  library(dplyr)
  
  #extract row containing enzymes
  e1=filter(relist, name == enzyme1)
  e2=filter(relist, name== enzyme2)
  
  
  #now extract restriction sites sequence 
  e1_sequence = as.character(tolower(e1$seq))
  e2_sequence = as.character (tolower(e2$seq))
  
 
  #locating restriction sites
  site1= str_locate_all(inputDNA,toupper(e1_sequence))
  site2=str_locate_all(inputDNA,toupper(e2_sequence))
  
  site1
})




output$site2 <-renderText({
  paste("Location of" ,input$enzyme2, "within input DNA sequence")}) 




output$search2 <-renderTable({
  
  #required inputs
  enzyme1 = input$enzyme1
  enzyme2= input$enzyme2
  inputDNA= input$inputDNA
  inputDNA=toupper(inputDNA)
  forward= input$forward
  reverse=input$reverse
  length= nchar(inputDNA)
  library(dplyr)
  library(dplyr)
  
  #extract row containing enzymes
  e1=filter(relist, name == enzyme1)
  e2=filter(relist, name== enzyme2)
  
  
  #now extract restriction sites sequence 
  e1_sequence = as.character(tolower(e1$seq))
  e2_sequence = as.character (tolower(e2$seq))
  
  
  #locating restriction sites
  site1= str_locate_all(inputDNA,toupper(e1_sequence))
  site2=str_locate_all(inputDNA,toupper(e2_sequence))
  
  site2
})

output$cut_list <-renderTable({
  #extract row containing enzymes
  enzyme1 = input$enzyme1
  enzyme2= input$enzyme2
  e1=filter(relist, name == enzyme1)
  e2=filter(relist, name== enzyme2)
  Enzymes= c(enzyme1,enzyme2)
  Site= c(as.character(e1$cut_site),as.character(e2$cut_site))
  cut_site_table= data.frame(Enzymes,Site)
  cut_site_table
})


}



shinyApp(ui,server)