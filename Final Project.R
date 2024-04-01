library(shiny)
library(deSolve)

# Step 1: Define the SEIR model function

SEIRmodel <- function(times, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * I * S
    dE <- beta * I * S - theta * E
    dI <- theta * E - gamma * I
    dR <- gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

server <- function(input, output) {
  # Step 2: Calculate SEIR model dynamics
  calculateSEIR <- function(beta, gamma, theta, S0, E0, I0, R0, tmax) {
    times <- seq(0, tmax, by = 1)
    parameters <- c(beta = beta, gamma = gamma, theta = theta)
    initialState <- c(S = S0, E = E0, I = I0, R = R0)
    solution <- ode(y = initialState, times = times, func = SEIRmodel, parms = parameters)
    return(solution)
  }
  
  # Step 3: Render the SEIR dynamics plot
  output$plot <- renderPlot({
    solution <- calculateSEIR(beta = input$transmission, gamma = input$recovery, theta = input$exposure, 
                              S0 = input$S0, E0 = input$E0, I0 = input$I0, R0 = input$R0, tmax = input$tmax)
    plot(solution[, "time"], solution[, "S"], type = "l", col = "black", lwd = 2, 
         xlab = "Time (days)", ylab = "Number of Individuals",
         main = "SEIR Model Dynamics")
    lines(solution[, "time"], solution[, "E"], col = "orange", lwd = 2)
    lines(solution[, "time"], solution[, "I"], col = "red", lwd = 2)
    lines(solution[, "time"], solution[, "R"], col = "blue", lwd = 2)
    legend("topright", legend = c("S", "E", "I", "R"), col = c("black", "orange", "red", "blue"), lwd = 2)
  })
  
  date_strings <- c("2021-10-16", "2021-10-23", "2021-10-30", "2021-10-31", 
                    "2021-11-07", "2021-11-14", "2021-11-21", "2021-11-28", 
                    "2021-12-05", "2021-12-12", "2021-12-19", "2021-12-26", 
                    "2022-01-02", "2022-01-09", "2022-01-16", "2022-01-23", 
                    "2022-01-30", "2022-02-06", "2022-02-13", "2022-02-20", 
                    "2022-02-27", "2022-03-06")
  
  
  # Step 4: Render the table of weekly cases
  output$weeklyCases <- renderTable({
    solution <- calculateSEIR(beta = input$transmission, gamma = input$recovery, theta = input$exposure, 
                              S0 = input$S0, E0 = input$E0, I0 = input$I0, R0 = input$R0, tmax = input$tmax)
    weekly_indices <- seq(1, length(solution[, "time"]), by = 7)
    I_weekly <- solution[weekly_indices, "I"]
    R_weekly <- solution[weekly_indices, "R"]
    
    new_infections_weekly <- c(I_weekly[1], diff(I_weekly)) + diff(c(0, R_weekly))
    
    formatted_dates <- format(as.Date(date_strings, format="%Y-%m-%d"), "%m-%d-%y")
    
    forecasted_data <- data.frame(Date = formatted_dates, New_Cases = head(new_infections_weekly, length(formatted_dates)))
    
    return(forecasted_data)
  })
  
  # Written by Arham
}

# Define UI
ui <- fluidPage(
  titlePanel("SEIR Model for Disease Spread with Dynamic Parameters"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("transmission", "Transmission Rate (beta):", min = 0, max = 1, value = 0.4234445346685245, step = 0.0001),
      sliderInput("recovery", "Recovery Rate (gamma):", min = 0, max = 1, value = 0.8511978960612641, step = 0.001),
      sliderInput("exposure", "Exposure Rate (theta):", min = 0, max = 0.01, value = 0.000435235234234, step = 0.00001),
      sliderInput("S0", "Initial Susceptible Individuals:", min = 100000, max = 10000000, value = 2040000, step = 10000),
      sliderInput("E0", "Initial Exposed Individuals:", min = 0, max = 1000, value = 590, step = 1),
      sliderInput("I0", "Initial Infected Individuals:", min = 1, max = 1000, value = 900, step = 1),
      sliderInput("R0", "Initial Recovered Individuals:", min = 0, max = 1000, value = 0, step = 10),
      sliderInput("tmax", "Time Span (days):", min = 30, max = 365, value = 154, step = 5)
    ),
    mainPanel(
      plotOutput("plot"),
      tableOutput("weeklyCases")
    )
  )
)

# Run the app
shinyApp(ui = ui, server = server)
