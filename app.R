############
# Function #
############
library(shiny)

Kappa.test.revised.fn<-function (x, y = NULL, null.value = 0, conf.level = 0.95) 
{
  DNAME <- deparse(substitute(x))
  METHOD <- paste0("Estimate Cohen's kappa statistics and test the null hypothesis ", 
                  "that kappa = ",null.value,".")
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (any(dim(x) < 2)) 
      stop("'x' must have at least 2 rows and columns")
    if (!is.numeric(x) || any(x < 0) || any(is.na(x))) 
      stop("all entries of 'x' must be nonnegative and finite")
  }
  else {
    if (is.null(y)) 
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y)) 
      stop("'x' and 'y' must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- factor(x[OK])
    y <- factor(y[OK])
    if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
      stop("'x' and 'y' must have at least 2 levels")
    x <- table(x, y)
  }
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr != nc) {
    stop("levels for 2 dimensions are different")
  }
  N <- sum(x)
  Po <- sum(diag(x))/N
  Pe <- sum(rowSums(x) * colSums(x)/N)/N
  kappa <- (Po - Pe)/(1 - Pe)
  JUDGEMENT <- c("No agreement", "Slight agreement", "Fair agreement", 
                 "Moderate agreement", "Substantial agreement", "Almost perfect agreement")
  judge <- JUDGEMENT[min(which(kappa < seq(0, 1, 0.2)))]
  seK0 <- sqrt(Pe/(N * (1 - Pe)))
  seK <- sqrt(Po * (1 - Po)/(N * (1 - Pe)^2))
  norm.pp <- qnorm(1 - (1 - conf.level)/2)
  Z <- (kappa - null.value)/seK0
  p.v <- 1 - pnorm(Z)
  kappaL <- kappa - norm.pp * seK
  kappaU <- kappa + norm.pp * seK
  CINT <- c(kappaL, kappaU)
  attr(CINT, "conf.level") <- conf.level
  RVAL <- list(statistic = c(Z = Z), estimate = kappa, conf.int = CINT, 
               p.value = p.v, method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  RVAL2 <- list(Result = RVAL, Judgement = judge, null.value = null.value)
  return(RVAL2)
}


###################
# ui file for app #
###################

ui <- fluidPage(
  
    titlePanel("Cohen's Kappa calculator"),
      sidebarLayout( fluid = TRUE, 
                 sidebarPanel(
                  numericInput("numcats", label = "Number of Categories", min = 2, max = 3, value = 2),
                  uiOutput("numbers"),
                  textInput("null.value", label = "Null value for Kappa", value = 0),
                  actionButton("actionButton", "Compute and Test Kappa!")
                 ),
                 
                 mainPanel(
                   verbatimTextOutput("text_out")
                 )
  )
)

#######################
# server file for app #
#######################

server <- function(input, output) {

  output$numbers<-renderUI({
    if(input$numcats==2){
      fluidRow(
        splitLayout(
          column(6,
                 numericInput("m1_cat1", "(1, 1)", value = NA),
                 numericInput("m1_cat2", "(1, 2)", value = NA)
          ),
          column(6,
                 numericInput("m2_cat1", "(2, 1)", value = NA),
                 numericInput("m2_cat2", "(2, 2)", value = NA)
          )
        )
      )
    }else{
      fluidRow(
        splitLayout(
          column(6,
                 numericInput("m1_cat1", "(1, 1)", value = NA),
                 numericInput("m1_cat2", "(1, 2)", value = NA),
                 numericInput("m1_cat3", "(1, 3)", value = NA)
          ),
          column(6,
                 numericInput("m2_cat1", "(2, 1)", value = NA),
                 numericInput("m2_cat2", "(2, 2)", value = NA),
                 numericInput("m2_cat3", "(2, 3)", value = NA)
          ),
          column(6,
                 numericInput("m3_cat1", "(3, 1)", value = NA),
                 numericInput("m3_cat2", "(3, 2)", value = NA),
                 numericInput("m3_cat3", "(3, 3)", value = NA)
          )
        )
      )
    }
  })
  
  results <- eventReactive(input$actionButton, {
  if(input$numcats==2){
    data<-matrix(c(input$m1_cat1, input$m1_cat2, 
                   input$m2_cat1, input$m2_cat2), nrow = 2, byrow = F)
  }else{
        data<-matrix(c(input$m1_cat1, input$m1_cat2, input$m1_cat3,
                       input$m2_cat1, input$m2_cat2, input$m2_cat3,
                       input$m3_cat1, input$m3_cat2, input$m3_cat3), nrow = 3, byrow = F)
      }
    null.value<-as.numeric(input$null.value)
    
    kappa_out<-Kappa.test.revised.fn(x = data, null.value = null.value)
    kappa_out
    })
  
  output$text_out <- renderPrint({
    results()
  })
  
}

#######################
# Run the application #
#######################

shinyApp(ui = ui, server = server)

#######
# END #
#######
