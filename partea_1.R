#a
frepcomgen <- function(m, n) {
  # Generați valorile lui X și Y
  x_values <- sort(sample(-(m^2):m^2, m, replace = FALSE))
  y_values <- sort(sample(-(n^2):n^2, n, replace = FALSE))
  

  
  pij_values <- matrix(runif(m * n), nrow = m, ncol = n)
  pij_values <- pij_values / sum(pij_values)
  
  # Generați valori suficiente pentru pi, qj și 𝜋ij
  
  qj_values <- colSums(pij_values)
  pi_values <- rowSums(pij_values)
  
  min_values <- min(m - 1, n - 1)
  indices_to_remove <- sample(seq_len(m * n), m+n - min_values )
  pij_values_flattened <- as.vector(pij_values)
  
  # Asigurați-vă că nu există două valori NA pe aceeași linie și coloană
  while (any(duplicated(indices_to_remove))) {
    indices_to_remove <- sample(seq_len(m * n), m+n - min_values )
  }
  pij_values_flattened[indices_to_remove] <- NA
  pij_values <- matrix(pij_values_flattened, nrow = m, ncol = n)
  
  
  matrice <- matrix(nrow = m+2, ncol = n+2)
  for (i in 1:m) {
    matrice[i + 1, 1] <- x_values[i]
  }
  
  for (j in 1:n) {
    matrice[1, j + 1] <- y_values[j]
  }
  
  for (i in 1:m) {
    for (j in 1:n) {
      matrice[i + 1, j + 1] <- pij_values[i, j]
    }
  }
  matrice[2:(m + 1), n + 2] <- pi_values
  
  # Add qj_values to the last row
  matrice[m + 2, 2:(n + 1)] <- qj_values
  
  matrice[m+2, n+2] <- 1
  
  matrice[m + 2, sample(2:(n + 1), min_values/2)] <- NA
  matrice[sample(2:(m + 1), min_values/2), n + 2] <- NA
  
  matrice[1, n+2] = 0
  matrice[m+2, 1] = 0
  
  
  return(matrice)
}

# Testați funcția frepcomgen pentru m = 2 și n = 3
repartitia_comuna_incompleta <- frepcomgen(3, 4)
print(repartitia_comuna_incompleta)


#b
fcomplrepcom <- function(matrice) {
  m <- nrow(matrice)
  n <- ncol(matrice)
  # Iterate until no more NA values exist
  while (sum(is.na(matrice)) != 1) {
    
    
    for (i in 2:m) {
      if (sum(is.na(matrice[i, 2:n])) == 1) {
        missing_col_index <- which(is.na(matrice[i, 2:n]))
        if (missing_col_index == n) {
          matrice[i, n] <- 1 - sum(matrice[i, 2:(n-1)], na.rm = TRUE)
        } else {
          matrice[i, missing_col_index+1] <- matrice[i, n ] - sum(matrice[i, 2:(n-1)], na.rm = TRUE)
        }
      }
    }
    
    for (j in 2:n) {
      if (sum(is.na(matrice[2:(m), j])) == 1) {
        missing_row_index <- which(is.na(matrice[2:m, j]))
        if (missing_row_index == m) {
          matrice[m, j] <- 1 - sum(matrice[2:(m-1), j], na.rm = TRUE)
        } else {
          matrice[missing_row_index+1, j] <- matrice[m, j ] - sum(matrice[2:(m-1), j], na.rm = TRUE)
        }
      }
    }
    
  }
  

  return(matrice)
}

repartitia_comuna <- fcomplrepcom(repartitia_comuna_incompleta)
print(repartitia_comuna)


#c
frepcomarginal <- function(matrice) {
  m <- nrow(matrice) - 2 # The number of X values
  n <- ncol(matrice) - 2 # The number of Y values

  marginal_x <- matrix(nrow = 2, ncol = m)
  marginal_y <- matrix(nrow = 2, ncol = n)

  marginal_x[1, ] <- matrice[2:(m + 1), 1]
  marginal_y[1, ] <- matrice[1, 2:(n + 1)]
  
  marginal_x[2, ] <- apply(matrice[2:(m + 1), 2:(n + 1)], 1, sum, na.rm = TRUE)
  marginal_y[2, ] <- apply(matrice[2:(m + 1), 2:(n + 1)], 2, sum, na.rm = TRUE)

  list(marginal_x = marginal_x, marginal_y = marginal_y)
}

repartitii_marginale <- frepcomarginal(repartitia_comuna)

marginal_x <- repartitii_marginale$marginal_x
marginal_y <- repartitii_marginale$marginal_y

print(marginal_x)
print(marginal_y)



#d
media <- function(matrix) {
  # Extragem valorile și probabilitățile
  values <- matrix[1, ]
  probabilities <- matrix[2, ]
  
  mean_value <- sum(values * probabilities)
  
  return(mean_value)
}

media2 <- function(matrix) {
  # Extragem valorile și probabilitățile
  values <- matrix[1, ]
  probabilities <- matrix[2, ]
  
  mean_value <- sum(values * (probabilities^2))
  
  return(mean_value)
}

media_produs <- function(matrix) {
  # Extragem valorile și probabilitățile pentru X și Y
  values_Y <- matrix[1, 2:(ncol(matrix)-1)]
  values_X <- matrix[2:(nrow(matrix)-1), 1]
  probabilities <- matrix[2:(nrow(matrix)-1), 2:(ncol(matrix)-1)]
  
  # Calculăm E[X * Y]
  expected_product <- sum(outer(values_X, values_Y, `*`) * probabilities)
  
  return(expected_product)
}


fpropcov <- function(a, b, c, d, marginal_x, marginal_y, repartitia_comuna) {
  
  E_X <- media(marginal_x)
  E_Y <- media(marginal_y)
  var_X <- media2(marginal_x) - (E_X^2)
  var_Y <- media2(marginal_y) - (E_Y^2)
  E_XY <- media_produs(repartitia_comuna)
  cov_XY <- E_XY - E_X*E_Y

  # Calculăm covarianța dintre Z și T folosind proprietățile covarianței
  cov_ZT <- a * c * cov_XY + b * d * cov_XY + a * d * var_Y + c * b * var_X
  
  return(cov_ZT)
}

# Exemplu de utilizare:
a <- 2
b <- 3
c <- 1
d <- 4

cov_ZT <- fpropcov(a, b, c, d, marginal_x, marginal_y, repartitia_comuna)
print(cov_ZT)


#e
fpcond <- function(matrice) {
  m <- nrow(matrice) - 2  # Subtracting the extra rows for totals
  n <- ncol(matrice) - 2  # Subtracting the extra columns for totals
  
  # P(X|Y) matrix initialization
  p_x_given_y <- matrix(nrow = m, ncol = n)
  # P(Y|X) matrix initialization
  p_y_given_x <- matrix(nrow = n, ncol = m)
  
  # Extract the pij values from the input matrix
  pij_values <- matrice[2:(m+1), 2:(n+1)]
  
  # Extract pi and qj values
  pi_values <- matrice[2:(m+1), n+2]
  qj_values <- matrice[m+2, 2:(n+1)]
  
  # Calculate P(X|Y) = P(X,Y) / P(Y)
  for (i in 1:m) {
    for (j in 1:n) {
      p_x_given_y[i, j] <- pij_values[i, j] / qj_values[j]
    }
  }
  
  # Calculate P(Y|X) = P(X,Y) / P(X)
  for (j in 1:n) {
    for (i in 1:m) {
      p_y_given_x[j, i] <- pij_values[i, j] / pi_values[i]
    }
  }
  
  # Combine the results in a list to return both conditional distributions
  conditional_distributions <- list(P_X_given_Y = p_x_given_y, P_Y_given_X = p_y_given_x)
  
  return(conditional_distributions)
}

conditional_probabilities <- fpcond(repartitia_comuna)
print(conditional_probabilities)


#f

# Functia pentru a extrage probabilitatea asociata valorilor lui x si y
fExtrageProbabilitateValori <- function(matrice, valoare_x, valoare_y) {
  # Gasirea indexului pentru valoarea x
  index_x <- which(matrice[2:(nrow(matrice) - 1), 1] == valoare_x)
  # Gasirea indexului pentru valoarea y
  index_y <- which(matrice[1, 2:(ncol(matrice) - 1)] == valoare_y)
  
  if (length(index_x) == 0 || length(index_y) == 0) {
    return (0)
  }
  
  # Extragerea probabilitatii
  probabilitate <- matrice[index_x + 1, index_y + 1]
  
  return(probabilitate)
}

# Exemplu de utilizare:
valoare_x <- 1 # Valoarea lui x din repartitia comuna
valoare_y <- 9 # Valoarea lui y din repartitia comuna
probabilitate <- fExtrageProbabilitateValori(repartitia_comuna, valoare_x, valoare_y)
print(paste("Probabilitatea asociata valorilor x =", valoare_x, "si y =", valoare_y, "este:", probabilitate))




#g

calculeaza_covarianta <- function(repartitia_comuna) {
  
  repartitii_marginale <- frepcomarginal(repartitia_comuna)
  
  marginal_x <- repartitii_marginale$marginal_x
  marginal_y <- repartitii_marginale$marginal_y
  
  E_X <- media(marginal_x)
  E_Y <- media(marginal_y)
  E_XY <- media_produs(repartitia_comuna)
  cov_XY <- E_XY - E_X*E_Y
  
  # Calculăm covarianța Cov(5X+9,-3Y-2)
  cov <- 5 * (-3) * cov_XY
  
  return(cov)
}

cov <- calculeaza_covarianta(repartitia_comuna)
print(cov)

probabilitatea_conditionata <- function(matrice, x_range, y_condition) {
  m <- nrow(matrice)
  n <- ncol(matrice) 
  # Extract the relevant columns and rows
  x_values <- matrice[2:(m - 1), 1]
  y_values <- matrice[1, 2:(n - 1)]
  pij_values <- matrice[2:m, 2:n]
  
  # Find indices of relevant values
  x_indices <- which(x_values > x_range[1] & x_values < x_range[2])
  y_condition_index <- which(y_values > y_condition)
  

  # Calculate the conditional probability
  numerator <- sum(pij_values[x_indices, y_condition_index])
  denominator <- sum(pij_values[m-1, y_condition_index])
  

  conditional_probability <- numerator / denominator
  return(conditional_probability)
}

# Set the range for X and the condition for Y
x_range <- c(0, 0.8)
y_condition <- 0.3

# Call the function with your joint distribution
conditional_probability <- probabilitatea_conditionata(repartitia_comuna, x_range, y_condition)
print(conditional_probability)



probabilitatea_interval <- function(matrice, x_great, y_lower) {
  m <- nrow(matrice)
  n <- ncol(matrice) 
  # Extract the relevant columns and rows
  x_values <- matrice[2:(m - 1), 1]
  y_values <- matrice[1, 2:(n - 1)]
  pij_values <- matrice[2:m, 2:n]
  
  # Find indices of relevant values
  x_indices <- which(x_values > x_great)
  y_condition_index <- which(y_values < y_lower)
  
  
  # Calculate the conditional probability
  probabilitatea <- sum(pij_values[x_indices, y_condition_index])

  return(probabilitatea)
}

# Set the range for X and the condition for Y
x_great <- 0.2
y_lower <- 1.7

# Call the function with your joint distribution
probabilitatea <- probabilitatea_interval(repartitia_comuna, x_great, y_lower)
print(probabilitatea)



#h
#1

# Functia pentru verificarea independentei variabilelor X si Y
verificaIndependenta <- function(matrice) {
  m <- nrow(matrice) - 2
  n <- ncol(matrice) - 2
  
  pi_values <- matrice[2:(m + 1), n + 2] # Probabilitatile marginale pentru X
  qj_values <- matrice[m + 2, 2:(n + 1)] # Probabilitatile marginale pentru Y
  
  for (i in 1:m) {
    for (j in 1:n) {
      prob_comuna <- matrice[i + 1, j + 1] # P(X=x, Y=y)
      prob_produs <- pi_values[i] * qj_values[j] # P(X=x) * P(Y=y)
      
      if (prob_comuna != prob_produs) {
        cat(sprintf("Variabilele X si Y nu sunt independente. Gasit la x=%d, y=%d\n", i, j))
        return(0)
      }
    }
  }
  
  cat("Variabilele X si Y sunt independente.\n")
  return(TRUE)
}

# Apelarea functiei cu matricea repartitiei comune ca argument
 verificaIndependenta(repartitia_comuna)


#2
calculeaza_corelatia <- function(repartitia_comuna) {
  
  repartitii_marginale <- frepcomarginal(repartitia_comuna)
  
  marginal_x <- repartitii_marginale$marginal_x
  marginal_y <- repartitii_marginale$marginal_y
  
  E_X <- media(marginal_x)
  E_Y <- media(marginal_y)
  E_XY <- media_produs(repartitia_comuna)
  
  var_X <- media2(marginal_x) - (E_X^2)
  var_Y <- media2(marginal_y) - (E_Y^2)
  
  cov_XY <- E_XY - E_X*E_Y
  
  # Calculăm covarianța Cov(5X+9,-3Y-2)
  corelatia_XY <- cov_XY/sqrt(var_X*var_Y)
  
  return(corelatia_XY)
}

corelatia_XY <- calculeaza_corelatia(repartitia_comuna)
print(corelatia_XY)

print(repartitia_comuna)

#i
m <- nrow(repartitia_comuna)
n <- ncol(repartitia_comuna) 
# Extract the relevant columns and rows
pij_values <- repartitia_comuna[3:m-1, 3:n-1]
print(pij_values)
m <- m-2
n <- n-2

Z <- c(1, 2, 3)

array_tridimensional <- array(0, dim = c(m, n, length(Z)))

# Distribuie probabilitatile din matrice pentru valorile lui Z
for (i in 1:m) {
  for (j in 1:n) {
    # Genereaza 3 probabilitati care insumate sunt egale cu probabilitatea din matrice
    probabilitati_Z <- runif(3)
    probabilitati_Z <- probabilitati_Z / sum(probabilitati_Z) * pij_values[i, j]
    
    # Atribuie valorile lui Z in functie de probabilitatile generate
    array_tridimensional[i, j, ] <- probabilitati_Z
    
    # Distribuie probabilitatile pentru fiecare valoare a lui Z
    # for (k in 1:length(Z)) {
    #   array_tridimensional[i, j, k] <- array_tridimensional[i, j, k] * probabilitati_Z[k]
    # }
  }
}

# Afiseaza array-ul tridimensional
print(array_tridimensional)
print(repartitia_comuna)



Z <- c(1, 2, 3)
x_values <- repartitia_comuna[2:(m - 1), 1]
y_values <- repartitia_comuna[1, 2:(n - 1)]
x_values
y_values



# Instalează și încarcă pachetul tidyverse, dacă nu este deja instalat
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)



# Crearea tabelului cu 3 coloane
data_table1 <- expand.grid(X = x_values, Y = y_values) %>%
  mutate(Z = as.vector(t(array_tridimensional[,,1])))

# Afișarea tabelului
print(data_table1)

data_table2 <- expand.grid(X = x_values, Y = y_values) %>%
  mutate(Z = as.vector(t(array_tridimensional[,,2])))

data_table3 <- expand.grid(X = x_values, Y = y_values) %>%
  mutate(Z = as.vector(t(array_tridimensional[,,3])))


plot3d(x = data_table1$X, y = data_table1$Y, z = data_table1$Z, col = "red", type = 's', radius = 0.2,
       xlab = "X", ylab = "Y", zlab = "Z")

# Add the second dataset to the plot
plot3d(x = data_table2$X, y = data_table2$Y, z = data_table2$Z, col = "blue", type = 's', radius = 0.2, add = TRUE)

# Add the third dataset to the plot
plot3d(x = data_table3$X, y = data_table3$Y, z = data_table3$Z, col = "yellow", type = 's', radius = 0.2, add = TRUE)

legend3d("topright", legend = c("Valorile lui Z:", "Rosu = 1", "Albastru = 2", "Galben = 3"), col = c("white", "red", "blue", "yellow"), pch = 16)