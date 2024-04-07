library(plotly)
library(shiny)
library(animate)
library(mvtnorm)
library(ggplot2)

calculate_double_integral <- function(f_expr, x_min, x_max, y_min, y_max) {
  f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  integral_xy <- try(integrate(Vectorize(function(x) {
    integrate(Vectorize(function(y) f(x, y)), y_min, y_max)$value
  }), x_min, x_max)$value, silent = TRUE)
  
  integral_yx <- try(integrate(Vectorize(function(y) {
    integrate(Vectorize(function(x) f(x, y)), x_min, x_max)$value
  }), y_min, y_max)$value, silent = TRUE)
  
  if (inherits(integral_xy, "try-error") || inherits(integral_yx, "try-error")) {
    return("Eroare la calculul integralei. Verificati daca functia este integrabila.")
  }
  
  if (abs(integral_xy - integral_yx) < 1e-5) {
    return(integral_xy)
  } else {
    return("Integrala nu poate fi calculata folosind Teorema lui Fubini.")
  }
}

is_probability_density <- function(f_expr, x_range, y_range) {
  f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  is_non_negative <- all(outer(x_range, y_range, Vectorize(f)) >= 0)
  integral_value <- calculate_double_integral(f_expr, min(x_range), max(x_range), min(y_range), max(y_range))
  is_integral_one <- abs(integral_value - 1) < 1e-5
  return(is_non_negative && is_integral_one)
}


#d
generate_continuous_random_variable <- function(f_expr, range_min, range_max, dimension) {
  f <- eval(parse(text = paste("function(", ifelse(dimension == 1, "x", "x, y"), ") {", f_expr, "}")))
  
  if (dimension == 1) {
    # Generare pentru variabila unidimensionala
    random_variable <- integrate(Vectorize(function(x) f(x)), range_min, range_max)$value
  } else {
    # Generare pentru variabila bidimensionala
    random_variable <- rcont2(n = 1, f = f, lower = c(range_min, range_min), upper = c(range_max, range_max))
  }
  
  return(random_variable)
}

# Functie pentru generarea variabilelor aleatoare continue
generate_continuous_rv <- function(pdf, n, dim = 1, lower = 0, upper = 1) {
  if (dim == 1) {
    # Pentru variabile aleatoare unidimensionale
    sample <- numeric(n)
    for (i in 1:n) {
      repeat {
        x <- runif(1, lower, upper)
        y <- runif(1, 0, max(pdf(x), na.rm = TRUE))
        if (y < pdf(x)) {
          sample[i] <- x
          break
        }
      }
    }
    return(sample)
  } else if (dim == 2) {
    # Pentru variabile aleatoare bidimensionale
    sample <- matrix(numeric(2 * n), nrow = n)
    for (i in 1:n) {
      repeat {
        x <- runif(1, lower[1], upper[1])
        y <- runif(1, lower[2], upper[2])
        z <- runif(1, 0, max(pdf(x, y), na.rm = TRUE))
        if (z < pdf(x, y)) {
          sample[i, ] <- c(x, y)
          break
        }
      }
    }
    return(sample)
  } else {
    stop("Dimensiunea specificata este incorecta. Utilizati 1 pentru unidimensional sau 2 pentru bidimensional.")
  }
}


#g
calculate_mean <- function(f_expr, x_range, y_range) {
  f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  mean_x <- integrate(Vectorize(function(x) x * integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value), min(x_range), max(x_range))$value
  mean_y <- integrate(Vectorize(function(y) y * integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value), min(y_range), max(y_range))$value
  mean <- c(mean_x, mean_y)
  return(mean)
}

calculate_variance <- function(f_expr, x_range, y_range) {
  f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  variance_x <- integrate(Vectorize(function(x) (x - calculate_mean(f_expr, x_range, y_range)[1])^2 * integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value), min(x_range), max(x_range))$value
  variance_y <- integrate(Vectorize(function(y) (y - calculate_mean(f_expr, x_range, y_range)[2])^2 * integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value), min(y_range), max(y_range))$value
  variance <- c(variance_x, variance_y)
  return(variance)
}

calculate_moments <- function(f_expr, x_range, y_range, order) {
  f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  moments_x <- sapply(1:order, function(k) integrate(Vectorize(function(x) x^k * integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value), min(x_range), max(x_range))$value)
  moments_y <- sapply(1:order, function(k) integrate(Vectorize(function(y) y^k * integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value), min(y_range), max(y_range))$value)
  moments <- list(moments_x = moments_x, moments_y = moments_y)
  return(moments)
}

calculate_centered_moments <- function(f_expr, x_range, y_range, order) {
  mean <- calculate_mean(f_expr, x_range, y_range)
  f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  centered_moments_x <- sapply(1:order, function(k) integrate(Vectorize(function(x) (x - mean[1])^k * integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value), min(x_range), max(x_range))$value)
  centered_moments_y <- sapply(1:order, function(k) integrate(Vectorize(function(y) (y - mean[2])^k * integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value), min(y_range), max(y_range))$value)
  centered_moments <- list(centered_moments_x = centered_moments_x, centered_moments_y = centered_moments_y)
  return(centered_moments)
}


ui <- fluidPage(
  titlePanel("Interpretarea geometrica a integralei duble si calculul integralei"),
  sidebarLayout(
    sidebarPanel(
      textInput("function_input", "Introduceti expresia functiei (ex: x^2 + y^2):"),
      numericInput("x_min_input", "Valoare minima pentru x:", 0),
      numericInput("x_max_input", "Valoare maxima pentru x:", 1),
      numericInput("y_min_input", "Valoare minima pentru y:", 0),
      numericInput("y_max_input", "Valoare maxima pentru y:", 1),
      actionButton("plot_button", "Genereaza Graficul si Calculeaza Integrala"),
      actionButton("density_button", "Calculeaza Densitati"),
      actionButton("mean_button", "Calculeaza media, dispersia"),
      radioButtons("dimension_radio", "Dimensiune variabila aleatoare:",
                   choices = c("Unidimensionala" = 1, "Bidimensionala" = 2),
                   selected = 1),
      
      conditionalPanel(
        condition = "input.dimension_radio == 1",
        numericInput("samples_uni_input", "Numar de esantioane (unidimensional):", 1000)
      ),
      
      conditionalPanel(
        condition = "input.dimension_radio == 2",
        numericInput("samples_bi_input", "Numar de esantioane (bidimensional):", 1000)
      ),
      actionButton("generate_samples_button", "Genereaza variabila aleatoare"),
      actionButton("generate_button", "Generarea grafurilor"),
      textInput("function_g_input", "Introduceti expresia functiei g(x):"),
      actionButton("calculate_statistics_button", "Calculeaza Medie si Dispersie"),
      verbatimTextOutput("mean_result_g"),
      verbatimTextOutput("variance_result_g"),
      actionButton("calculate_probabilities_button", "Calculeaza Probabilitati"),
      actionButton("calculate_cov_corr_button", "Calculeaza Covarianta si Corelatia"),
      verbatimTextOutput("covariance_result"),
      verbatimTextOutput("correlation_result")
    ),
    mainPanel(
      plotlyOutput("plot"),
      textOutput("integral_result"),
      textOutput("density_check"),
      textOutput("marginal_density_x_result"),
      textOutput("marginal_density_y_result"),
      textOutput("conditional_density_y_given_x_result"),
      textOutput("conditional_density_x_given_y_result"),
      textOutput("mean_result"),
      textOutput("variance_result"),
      textOutput("moments_result"),
      textOutput("centered_moments_result"),
      textOutput("random_variable_output"),
      plotOutput("plot_samples"),
      textOutput("error_message"),
      plotOutput("density_plot"),
      plotOutput("cdf_plot"),
      textOutput("marginal_prob_X_result"),
      textOutput("marginal_prob_Y_result"),
      textOutput("conditional_prob_Y_given_X_result"),
      textOutput("conditional_prob_X_given_Y_result")
      # textOutput("mean_result_g"),
      # textOutput("variance_result_g"),
      # textOutput("covariance_result"),
      # textOutput("correlation_result")
      
    )
  )
)

server <- function(input, output) {
  calculate_double_integral <- function(f_expr, x_min, x_max, y_min, y_max) {
    f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
    integral_xy <- try(integrate(Vectorize(function(x) {
      integrate(Vectorize(function(y) f(x, y)), y_min, y_max)$value
    }), x_min, x_max)$value, silent = TRUE)
    
    integral_yx <- try(integrate(Vectorize(function(y) {
      integrate(Vectorize(function(x) f(x, y)), x_min, x_max)$value
    }), y_min, y_max)$value, silent = TRUE)
    
    if (inherits(integral_xy, "try-error") || inherits(integral_yx, "try-error")) {
      return("Eroare la calculul integralei. Verificati daca functia este integrabila.")
    }
    
    if (abs(integral_xy - integral_yx) < 1e-5) {
      return(integral_xy)
    } else {
      return("Integrala nu poate fi calculata folosind Teorema lui Fubini.")
    }
  }
  
  is_probability_density <- function(f_expr, x_range, y_range) {
    f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
    is_non_negative <- all(outer(x_range, y_range, Vectorize(f)) >= 0)
    integral_value <- calculate_double_integral(f_expr, min(x_range), max(x_range), min(y_range), max(y_range))
    is_integral_one <- abs(integral_value - 1) < 1e-5
    return(is_non_negative && is_integral_one)
  }
  
  calculate_marginal_densities <- function(f_expr, x_range, y_range) {
    f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
    f_x <- function(x) integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value
    f_y <- function(y) integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value
    list(f_x = f_x, f_y = f_y)
  }
  
  calculate_conditional_densities <- function(f_expr, x_range, y_range) {
    f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
    f_y_given_x <- function(y, x) f(x, y) / integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value
    f_x_given_y <- function(x, y) f(x, y) / integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value
    list(f_y_given_x = f_y_given_x, f_x_given_y = f_x_given_y)
  }
  
  output$plot <- renderPlotly({
    req(input$plot_button)
    f_expr <- input$function_input
    f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
    if (!is.function(f)) {
      stop("Expresia introdusa nu este o functie valida.")
    }
    
    x <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    
    z <- outer(x, y, Vectorize(f))
    
    plot_ly(x = ~x, y = ~y, z = ~z, type = "surface")
  })
  
  output$integral_result <- renderText({
    req(input$plot_button)
    f_expr <- input$function_input
    integral_result <- calculate_double_integral(f_expr, input$x_min_input, input$x_max_input, input$y_min_input, input$y_max_input)
    if (is.numeric(integral_result)) {
      paste("Valoarea integralei duble este:", integral_result)
    } else {
      integral_result
    }
  })
  
  output$density_check <- renderText({
    req(input$plot_button)
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    if (is_probability_density(f_expr, x_range, y_range)) {
      "Functia introdusa este o densitate de probabilitate valida."
    } else {
      "Functia introdusa NU este o densitate de probabilitate valida."
    }
  })
  
  # output$random_variable_output <- renderPrint({
  #   req(input$density_button)
  #   f_expr <- input$function_input
  #   x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
  #   y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
  #   
  #   if (is_probability_density(f_expr, x_range, y_range)) {
  #     "Generare variabila aleatoare..."
  #     f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
  #     random_variable <- r2dtable(n = 1, prob = outer(x_range, y_range, Vectorize(f)))
  #     random_variable
  #   } else {
  #     "Functia introdusa NU este o densitate de probabilitate valida. Nu se poate genera variabila aleatoare."
  #   }
  # })
  
  
  output$density_plot <- renderPlot({
    req(input$generate_button)
    f_expr <- input$function_input
    range_min <- input$range_min_input
    range_max <- input$range_max_input
    dimension <- input$dimension_radio
    
    x_vals <- seq(range_min, range_max, length.out = 100)
    
    if (dimension == 1) {
      f <- eval(parse(text = paste("function(x) {", f_expr, "}")))
      density_vals <- sapply(x_vals, function(x) f(x))
    } else {
      f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
      density_vals <- sapply(x_vals, function(x) integrate(Vectorize(function(y) f(x, y)), range_min, range_max)$value)
    }
    
    ggplot() +
      geom_line(aes(x = x_vals, y = density_vals), color = "blue") +
      labs(title = "Densitatea de probabilitate", x = "Valoarea variabilei", y = "Densitatea")
  })
  
  output$cdf_plot <- renderPlot({
    req(input$generate_button)
    f_expr <- input$function_input
    range_min <- input$range_min_input
    range_max <- input$range_max_input
    dimension <- input$dimension_radio
    
    x_vals <- seq(range_min, range_max, length.out = 100)
    
    if (dimension == 1) {
      f <- eval(parse(text = paste("function(x) {", f_expr, "}")))
      cdf_vals <- sapply(x_vals, function(x) integrate(Vectorize(f), range_min, x)$value)
    } else {
      f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
      cdf_vals <- sapply(x_vals, function(x) integrate(Vectorize(function(y) f(x, y)), range_min, x)$value)
    }
    
    ggplot() +
      geom_line(aes(x = x_vals, y = cdf_vals), color = "red") +
      labs(title = "Functia de repartitie", x = "Valoarea variabilei", y = "Probabilitatea cumulata")
  })
  
  
  output$mean_result <- renderText({
    req(input$mean_button)
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    mean <- calculate_mean(f_expr, x_range, y_range)
    paste("Media pentru X:", mean[1], "Media pentru Y:", mean[2])
  })
  
  output$variance_result <- renderText({
    req(input$mean_button)
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    variance <- calculate_variance(f_expr, x_range, y_range)
    paste("Dispersia pentru X:", variance[1], "Dispersia pentru Y:", variance[2])
  })
  
  output$moments_result <- renderText({
    req(input$mean_button)
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    order <- 4  # Puteti schimba ordinul la care doriti sa calculati momentele
    moments <- calculate_moments(f_expr, x_range, y_range, order)
    paste("Momente pentru X:" , moments$moments_x, "Momente pentru Y:", moments$moments_y)
  })
  
  output$centered_moments_result <- renderText({
    req(input$mean_button)
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    order <- 4  # Puteti schimba ordinul la care doriti sa calculati momentele centrate
    centered_moments <- calculate_centered_moments(f_expr, x_range, y_range, order)
    paste("Momente centrate pentru X:", centered_moments$centered_moments_x, "Momente centrate pentru Y:", centered_moments$centered_moments_y)
  })
  
  output$random_variable_output <- renderPrint({
    req(input$generate_button)
    f_expr <- input$function_input
    range_min <- input$range_min_input
    range_max <- input$range_max_input
    dimension <- input$dimension_radio
    
    if (dimension == 1) {
      random_variable <- generate_continuous_random_variable(f_expr, range_min, range_max, dimension)
      paste("Valoarea variabilei aleatoare unidimensionale generate este:", random_variable)
    } else {
      random_variable <- generate_continuous_random_variable(f_expr, range_min, range_max, dimension)
      paste("Valoarea variabilei aleatoare bidimensionale generate este:", random_variable)
    }
  })
  
  observeEvent(input$density_button, {
    req(input$function_input)
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    
    marginal_densities <- calculate_marginal_densities(f_expr, x_range, y_range)
    conditional_densities <- calculate_conditional_densities(f_expr, x_range, y_range)
    
    output$marginal_density_x_result <- renderText({
      paste("Densitatea marginala pentru X:", marginal_densities$f_x(1))
    })
    
    output$marginal_density_y_result <- renderText({
      paste("Densitatea marginala pentru Y:", marginal_densities$f_y(1))
    })
    
    output$conditional_density_y_given_x_result <- renderText({
      paste("Densitatea conditionata pentru Y dat X:", conditional_densities$f_y_given_x(1, 0.5))
    })
    
    output$conditional_density_x_given_y_result <- renderText({
      paste("Densitatea conditionata pentru X dat Y:", conditional_densities$f_x_given_y(0.5, 1))
    })
  })
  
  observeEvent(input$generate_samples_button, {
    req(input$function_input)
    
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    
    if (is_probability_density(f_expr, x_range, y_range)) {
      if (input$dimension_radio == 1) {
        # Generare variabila aleatoare unidimensionala
        samples_uni <- generate_continuous_rv(eval(parse(text = paste("function(x) {", f_expr, "}"))),
                                              n = input$samples_uni_input, dim = 1,
                                              lower = input$x_min_input, upper = input$x_max_input)
        # Afisare unidimensionala
        output$plot_samples <- renderPlot({
          hist(samples_uni, breaks = 30, main = "Esantioane din distributia normala unidimensionala")
        })
      } else if (input$dimension_radio == 2) {
        # Generare variabila aleatoare bidimensionala
        samples_bi <- generate_continuous_rv(eval(parse(text = paste("function(x, y) {", f_expr, "}"))),
                                             n = input$samples_bi_input, dim = 2,
                                             lower = c(input$x_min_input, input$y_min_input),
                                             upper = c(input$x_max_input, input$y_max_input))
        # Afisare bidimensionala
        output$plot_samples <- renderPlot({
          plot(samples_bi[, 1], samples_bi[, 2], xlab = "X", ylab = "Y",
               main = "Esantioane din distributia normala bidimensionala")
        })
      }
    } else {
      # Afisare un mesaj de eroare daca functia nu este o densitate de probabilitate valida
      # Poti afisa mesajul intr-un element de tip textOutput in ui.R
      output$error_message <- renderText("Functia introdusa NU este o densitate de probabilitate valida.")
      output$plot_samples <- renderPlot(NULL)  # Se sterge graficul daca exista o eroare
    }
  })
  
  observeEvent(input$calculate_statistics_button, {
    req(input$function_g_input)
    
    g_expr <- input$function_g_input
    f_X_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    
    # Defineste functia g(x)
    g <- eval(parse(text = paste("function(x) {", g_expr, "}")))
    
    # Defineste functia de densitate de probabilitate f_X(x) pentru variabila X
    f_X <- eval(parse(text = paste("function(x) {", f_X_expr, "}")))
    
    # Calculeaza media E[g(X)]
    mean_g_X <- integrate(function(x) g(x) * f_X(x), lower = input$x_min_input, upper = input$x_max_input)$value
    
    # Calculeaza dispersia Var[g(X)]
    variance_g_X <- integrate(function(x) (g(x) - mean_g_X)^2 * f_X(x), lower = input$x_min_input, upper = input$x_max_input)$value
    
    # Afiseaza rezultatele
    output$mean_result_g <- renderText({
      paste("Media E[g(X)] este:", mean_g_X)
    })
      
    output$variance_result_g <- renderText({
      paste("Dispersia Var[g(X)] este:", variance_g_X)
    })
  })
  
  
  # Functie pentru calculul diferitelor tipuri de probabilitati
  calculate_probabilities <- function(f_expr, x_range, y_range) {
    f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
    
    # Verifica daca functia este o densitate de probabilitate valida
    if (!is_probability_density(f_expr, x_range, y_range)) {
      return("Functia introdusa NU este o densitate de probabilitate valida.")
    }
    
    # Calculeaza probabilitatea marginala pentru X
    marginal_prob_X <- function(x) integrate(Vectorize(function(y) f(x, y)), min(y_range), max(y_range))$value
    
    # Calculeaza probabilitatea marginala pentru Y
    marginal_prob_Y <- function(y) integrate(Vectorize(function(x) f(x, y)), min(x_range), max(x_range))$value
    
    # Calculeaza probabilitatea conditionata pentru Y dat X
    conditional_prob_Y_given_X <- function(y, x) f(x, y) / marginal_prob_X(x)
    
    # Calculeaza probabilitatea conditionata pentru X dat Y
    conditional_prob_X_given_Y <- function(x, y) f(x, y) / marginal_prob_Y(y)
    
    # Returneaza rezultatele sub forma de lista
    list(marginal_prob_X = marginal_prob_X,
         marginal_prob_Y = marginal_prob_Y,
         conditional_prob_Y_given_X = conditional_prob_Y_given_X,
         conditional_prob_X_given_Y = conditional_prob_X_given_Y)
  }
  
  # Adauga aceasta functie in interiorul server <- function(input, output) { ... }
  
  observeEvent(input$calculate_probabilities_button, {
    req(input$function_input)
    
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    
    probabilities <- calculate_probabilities(f_expr, x_range, y_range)
    
    # Afiseaza rezultatele  
    output$marginal_prob_X_result <- renderPrint({
      paste("Probabilitatea marginala pentru X:", probabilities$marginal_prob_X(1))
    })
    
    output$marginal_prob_Y_result <- renderPrint({
      paste("Probabilitatea marginala pentru Y:", probabilities$marginal_prob_Y(1))
    })
    
    output$conditional_prob_Y_given_X_result <- renderPrint({
      paste("Probabilitatea conditionata pentru Y dat X:", probabilities$conditional_prob_Y_given_X(1, 0.5))
    })
    
    output$conditional_prob_X_given_Y_result <- renderPrint({
      paste("Probabilitatea conditionata pentru X dat Y:", probabilities$conditional_prob_X_given_Y(0.5, 1))
    })
  })
  
  
  
  
  
  
  observeEvent(input$calculate_cov_corr_button, {
    req(input$function_input)
    
    f_expr <- input$function_input
    x_range <- seq(input$x_min_input, input$x_max_input, length.out = 50)
    y_range <- seq(input$y_min_input, input$y_max_input, length.out = 50)
    
    if (is_probability_density(f_expr, x_range, y_range)) {
      # Defineste functia de densitate de probabilitate f(x, y)
      f <- eval(parse(text = paste("function(x, y) {", f_expr, "}")))
      
      # Genereaza esantioane din f(x, y)
      samples_bi <- generate_continuous_rv(f, n = 1000, dim = 2,
                                           lower = c(input$x_min_input, input$y_min_input),
                                           upper = c(input$x_max_input, input$y_max_input))
      
      # Calculeaza Covarianta si Coeficientul de Corelatie
      cov_XY <- cov(samples_bi[, 1], samples_bi[, 2])
      cor_XY <- cor(samples_bi[, 1], samples_bi[, 2])
      
      # Afiseaza rezultatele
      output$covariance_result <- renderPrint({
        paste("Covarianta Cov[X, Y] este:", cov_XY)
      })
      
      output$correlation_result <- renderPrint({
        paste("Coeficientul de corelatie ρ[X, Y] este:", cor_XY)
      })
    } else {
      # Afiseaza un mesaj de eroare daca functia nu este o densitate de probabilitate valida
      output$covariance_result <- renderPrint("Functia introdusa NU este o densitate de probabilitate valida.")
      output$correlation_result <- renderPrint("")
    }
  })
  
  
}

shinyApp(ui, server)
