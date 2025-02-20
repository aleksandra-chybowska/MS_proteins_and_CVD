library(tidyverse)
library(data.table)

setwd("~/Projects/python/proteins/data/phenotypes/")
proteins <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")
# the input dataset should contain 15,818 observations of 440 variables
dim(proteins)
# exclude "id" col. 439 protein IDs left.
cols = colnames(proteins)[-1] 
# find all protein groups: 306 items
groups = cols[grepl('\\.', cols)] 

# if we remove protein groups from all protein names, 
# we will get 133 IDs of individual proteins
# 133 + 306 = 439
individual = setdiff(cols, groups) 

# Now focus on protein groups
# Split protein group names by dot
prot_names = unlist(strsplit(groups, '\\.'))
# Find unique IDs
prot_names = unique(prot_names)
# Finally, concatenate individual and protein group IDs into one vector
all = c(individual, prot_names)
# you should get 596 unique IDs
all = unique(all) 

# paste all individual proteins into one string that you can pass to
# Retrieve / ID mapping UniProt tool
all_str = paste(all, collapse=' ')

# We paste the generated protein ID list into Retrieve/ ID mapping tool 
# available here: https://www.uniprot.org/id-mapping, download the results
# and read them in
setwd("~/Projects/python/proteins/data/annotations/ids/")
annots = read.csv("idmapping_2025_01_20.tsv", sep='\t')
# Clean up protein names. This procedure leaves some deleted proteins in 
# but it is fine atm
annots$Protein.names = tstrsplit(annots$Protein.names, "\\(")[[1]]
annots$Protein.names = tstrsplit(annots$Protein.names, "\\[")[[1]]
annots$Protein.names = str_trim(annots$Protein.names)


# First, we will annotate individual proteins
unique_df = data.frame("id"=individual)
unique_df = merge(unique_df, annots[c("From", "Protein.names")], 
                  by.x="id", by.y="From")
# Final name
unique_df$Name = unique_df$Protein.names


# OK, now the harder part. Protein groups.
# Get IDs of all protein groups
groups_df = data.frame("id"=groups)

# Now I need to translate A.B.C into three rows, A, B, and C
groups_df = groups_df %>%
  mutate(components = str_split(id, "\\.")) %>% 
  unnest_longer(components)

groups_df$order = 1:nrow(groups_df)
# Merge with annots
groups_df = merge(groups_df, annots[c("From", "Protein.names", "Gene.Names")], 
                  by.x="components", by.y="From") %>%
            select(id, components, Protein.names, Gene.Names, order) %>%
            arrange(order)

# Divide into groups and mixed groups
# Group protein groups by ID. Ignore deleted records. Find same genes.
groups_groupped = groups_df %>%
  group_by(id) %>%
  filter(Protein.names != "deleted") %>%
  summarise(
    all_same_name = n_distinct(Protein.names) == 1,
    all_same_gene = {
        # Split each Gene.Names into individual words
        gene_lists <- strsplit(Gene.Names, " ")
        
        # Find the intersection of all gene lists within the group
        common_genes <- Reduce(intersect, gene_lists)
        
        # Check if there's at least one common gene across all rows
        length(common_genes) > 0
      }
    )

groups_df = merge(groups_df, groups_groupped, by="id", all.x=TRUE)


# Nice! Now, get groups of isoforms. Add a (G) to their name. 
groups_df$type = ifelse(groups_df$all_same_gene , "G", "PG")
# write.csv(groups_df, "Debugging_DS.csv")

# read a file with the current order of PGs, leave only PGs
order = read.csv("order.csv") %>%
  filter(!id %in% individual)
# wrong order, correct
order_pg = order %>%
  filter(Name %like% "PG") %>%
  arrange(as.numeric(gsub("PG","", Name)))
order_g = order %>%
  filter(!(Name %like% "PG"))

correct_order = rbind(order_g, order_pg)

groups_named = groups_df %>%
  group_by(id) %>%
  summarise(
    Gene = str_c(na.omit(Gene.Names), collapse = ", "),
    Proteins = str_c(na.omit(Protein.names), collapse = ", "),
    Type = max(type),
    Name = {
      names_list <- na.omit(Protein.names)
      if (length(names_list) > 1 && names_list[1] == "deleted") {
        names_list[2]
      } else {
        first(names_list)
      }
    }
  ) %>%
  arrange(match(id, correct_order$id))


PGs = groups_named %>%
  filter(Type == "PG") %>%
  mutate(
    Numbered_Type = {
        type_counts <- ave(Type, Type, FUN = seq_along)
        paste0("PG ", type_counts)
      }
    )

Gs = groups_named %>%
  filter(Type == "G") %>%
  mutate(
    Numbered_Type = {
        type_counts <- ave(Name, Name, FUN = seq_along)
        total_counts <- ave(Name, Name, FUN = length)
        ifelse(total_counts == 1, "G", paste0("G", type_counts))
      }
    )
  
# combine group annots into a single df
group_annots = rbind(Gs, PGs)

# select only name from unique 
short_annots = group_annots %>% 
  mutate( 
    Full_name = case_when(
      Type == "PG" ~ Numbered_Type, 
      TRUE ~ paste0(Name, " (", Numbered_Type, ")")
    )) %>%
  select(id, Full_name)

colnames(short_annots) = c("id", "Name")

short_annots = rbind(unique_df[c("id", "Name")], short_annots)

## test with Hannah
manual_hannah = read.csv("U:/GS_mass_spec_proteins/phenotypes/short_annots_final.csv")
manual_and_gen = merge(manual_hannah, short_annots, by="id")
manual_and_gen$match <- manual_and_gen$Name.x == manual_and_gen$Name.y

table(manual_and_gen$Name.x, manual_and_gen$Name.y)
# 
# gen = read.csv("Protein_groups_annots_debug_manual.csv")
# gen = merge(manual_and_gen, 
#             gen[c("Original_id", "Type_Ola", "Suggested_new_name_Ola")], 
#             by.x="id", by.y = "Original_id")
# 
# write.csv(gen, "first_annots_debug_hannah.csv")
# ## test with Josie
# manual = read.csv("Josie_manual_annots.csv") %>%
#   filter(!original_id %in% individual) 
# 
# manual_groupped = manual %>%
#   group_by(original_id) %>%
#   summarise(
#     current_name = first(current_name), 
#     uniprot_gene = str_c(na.omit(uniprot_gene), collapse = ", "),
#     uniprot_name = str_c(na.omit(uniprot_name), collapse = ", "),
#     current_name_ok = first(na.omit(current_name_ok)),
#     suggested_new_name = first(na.omit(Suggested_new_name))
#   ) %>%
#   ungroup()
# 
# manual_groupped <- manual_groupped %>%
#   arrange(match(original_id, order$id))
# 
# groups_joined = left_join(manual_groupped, groups_named, 
#                           by=c("original_id" = "id"))
# colnames(groups_joined) = c("Original_id",
#                             "Current_name",
#                             "Gene_Josie",
#                             "Name_Josie",
#                             "Current_name_ok_Josie",
#                             "Suggested_new_name_Josie",
#                             "Gene_Ola",
#                             "Name_Ola",
#                             "Type_Ola",
#                             "Suggested_new_name_Ola")
# write.csv(groups_joined, "Protein_groups_annots_debug.csv", row.names = F)
