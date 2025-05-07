library(httr)
library(jsonlite)
library(tidyverse)

# Define the GraphQL API endpoint

diseases <- list(
  "composite_CVD" = "EFO_0000319",
  "hf" = "EFO_0003144",
  "isch_stroke" = "HP_0002140", # "EFO_0000712",
  "mi" = "EFO_0000612",
  "tia" = "EFO_0003764",
  "chd" = "EFO_0001645",
  "death" = "EFO_0005056"
)

# Define the GraphQL queries

fetch_disease <- function(efoID, count) {
  query <- paste0("{
    disease(efoId: \"", efoID ,"\") {
      associatedTargets(page: {size: ", count,", index:0}) {
        rows {
          target {
            id
            approvedSymbol
            proteinIds {
              id
            } 
          }
        }
      }
    }
  }")
  return(query)
}

fetch_disease_count <- function(efoID) {
  query <- paste0("{
    disease(efoId: \"", efoID ,"\") {
      associatedTargets(page: {size: 10, index:0}) {
        count
        rows {
          target {
            id
          }
        }
      }
    }
  }")
  return(query)
}

fetch_protein <- function(uniprotID, count) {
  query <- paste0("{
              mapIds(queryTerms: [", uniprotID, "]){
                mappings{
                  hits{
                    id
                    name
                    object {
                      ... on Target {
                        id
                        associatedDiseases(page: {size: ", count, ", index:0}) {
                          count
                          rows {
                            disease {
                              id
                              name
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }")

  return(query)
}

fetch_protein_count <- function(uniprotID) {
  query <- paste0("{
              mapIds(queryTerms: [", uniprotID, "]){
                mappings{
                  hits{
                    id
                    name
                    object {
                      ... on Target {
                        id
                        associatedDiseases(page: {size: 10, index:0}) {
                          count
                          rows {
                            disease {
                              id
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }")
  
  return(query)
}

run_OT_query <- function(query) {
  # Send POST request to Open Targets API
  url <- "https://api.platform.opentargets.org/api/v4/graphql"
  response <- POST(url, body = list(query = query), encode = "json")
  
  # Parse the response
  if (status_code(response) == 200) {
    data <- content(response, as = "text", encoding = "UTF-8")
    parsed_data <- fromJSON(data, flatten = TRUE)
  } else {
    print(paste("Error:", status_code(response)))
  }
  
  return(parsed_data)
}

disease_count <- function(efoID) {
  query <- fetch_disease_count(efoID)
  res <- run_OT_query(query)
  
  return(res$data$disease$associatedTargets$count)
}

find_max_disease_count <- function(diseases) {
  # max = 0
  # for (disease in diseases){
  #   count = disease_count(disease)
  #   if (count > max) {
  #     max = count
  #   }
  # }
  # return(max)
  
  max_count <- max(sapply(diseases, disease_count), na.rm = TRUE)
  return(max_count)

}

get_disease <- function(efoID, max_count) {
  query <- fetch_disease(efoID, max_count)
  res <- run_OT_query(query)
  
  return(res$data$disease$associatedTargets$rows)
}

parse_disease <- function(efoID, max_count, disease_name) {
  disease <- get_disease(efoID, max_count)
  disease_df <- as_tibble(disease) 
  
  # Unnest the protein IDs column
  disease_tidy <- disease_df %>% unnest(target.proteinIds)
  colnames(disease_tidy) <- c("Ensembl_ID", "ApprovedSymbol", "UniProt_ID")
  disease_tidy$Disease <- disease_name
  return(disease_tidy)
}

# main  
setwd("~/Projects/python/proteins/data/disease/open_targets")

max_disease_count <- find_max_disease_count(diseases)
cat(paste0("Max count:  ", max_disease_count,"\n"))
disease_keys <- names(diseases)
disease_data <- list()

# Loop through each disease key, parse data, and save results
for (key in disease_keys) {
  cat(paste0("Processing: ", key,"\n"))
  disease_data[[key]] <- parse_disease(diseases[[key]], max_disease_count, key)
  
  # Save each dataframe as a CSV file
  write.csv(disease_data[[key]], paste0(key, "_ot_data.csv"), row.names = FALSE)
  # saveRDS(disease_data[[key]], paste0(key, "_ot_data.RDS"))
}

# Print confirmation
cat("Processing complete. Data saved as CSV files.\n")
