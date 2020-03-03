# --------------- #
#### Functions ####
# --------------- #

require(data.table)
require(reshape2)
require(MPTinR)
require(tidyverse)


#' Preparing Data for ReAL Model estimation
#'
#' This function allows you to prepare IAT data for the application of functions estimating ReAL model parameters.
#' @param data Dataframe containing all relevant columns (Subject, TaskSwitch, Trial, Compatibility, StimulusCat and Correct).
#' @param Subject Name of the column containing the subject id.
#' @param TaskSwitch Name of the column containing infomation about task switching.
#' @param Switch Level of TaskSwitch coding task switchung.
#' @param Repetition Level of TaskSwitch coding task repetition.
#' @param Trial Name of the column containing the trial number.
#' @param Compatibility Name of the column containing the compatibility of current key mapping.
#' @param Compatible Level of Compatibility coding a compatible mapping.
#' @param Incompatible Level of Compatibility coding an incompatible mapping.
#' @param StimulusCat Name of the column containing the category of the current stimulus.
#' @param Correct Name of the column containing accuracy of responses (0 = error, 1 = correct).
#' @param TargetCat1 Level of StimulusCat coding the first target category (i.e, flower)
#' @param TargetCat2 Level of StimulusCat coding the second target category (i.e, insect)
#' @param AttributeCat1 Level of StimulusCat coding the first attribute category (i.e, pleasant)
#' @param AttributeCat2 Level of StimulusCat coding the second attribute category (i.e, unpleasant)
#' @keywords IAT mpt ReAL
#' @export
#' @examples
#' prep_ReAL()
#'
prep_ReAL <- function(data,
                      Subject,
                      TaskSwitch = NULL,
                      Switch = 1,
                      Repetition = 0,
                      Trial = NULL,
                      Compatibility,
                      Compatible = 1,
                      Incompatible = 2,
                      StimulusCat,
                      Correct,
                      TargetCat1,
                      TargetCat2,
                      AttributeCat1,
                      AttributeCat2){

  if(class(data)[1] != "data.frame") stop('Data must be a data.frame!')
  if(all(c(Subject,Compatibility,StimulusCat,Correct)  %in% names(data)) == FALSE  ) stop('Check column names!')
  if(setequal(c(TargetCat1,TargetCat2,AttributeCat1,AttributeCat2), unique(data[[StimulusCat]])) == FALSE  ) stop('Check category names!')
  if(is.null(TaskSwitch) & is.null(Trial)) stop('Specify either TaskSwitch or Trial!')


  tmp <- data.table(Subject = data[[Subject]],
                    Compatibility = factor(data[[Compatibility]], levels = c(1,2), labels = c("Compatible","Incompatible")),
                    Accuracy = data[[Correct]])

  tmp$StimulusCat[data[[StimulusCat]] == TargetCat1] <- "TargetCat1"
  tmp$StimulusCat[data[[StimulusCat]] == TargetCat2] <- "TargetCat2"

  tmp$StimulusCat[data[[StimulusCat]] == AttributeCat1] <- "AttributeCat1"
  tmp$StimulusCat[data[[StimulusCat]] == AttributeCat2] <- "AttributeCat2"

  # Task Switch
  if(is.null(TaskSwitch)){

    tmp$Task <- ifelse(tmp$StimulusCat == "AttributeCat1" | tmp$StimulusCat == "AttributeCat2", "Attribute", "Target")

    tmp$TaskSwitch <- numeric()

    # vectorization
    tmp$TaskSwitch <- ifelse(tmp$Task == c(NA,tmp$Task)[-(length(tmp$Task)+1)], "TR", "TS")
    tmp$TaskSwitch[data[[Trial]] < c(NA,data[[Trial]])[-(length(data[[Trial]])+1)]] <- NA

  } else{

    tmp$TaskSwitch <- data[[TaskSwitch]]
    tmp$TaskSwitch[data[[TaskSwitch]] == Switch] <- "TS"
    tmp$TaskSwitch[data[[TaskSwitch]] == Repetition] <- "TR"

  }

  ReAL_data <- tmp[is.na(TaskSwitch) == FALSE,
                   .(correct = sum(Accuracy == 1),
                     incorrect = sum(Accuracy == 0)),
                   by = c("Subject", "StimulusCat", "Compatibility", "TaskSwitch")]

  # -> WIDE FORMAT
  ReAL_data <- dcast.data.table(ReAL_data, Subject ~ StimulusCat + TaskSwitch + Compatibility, value.var = c("incorrect","correct"), sep = " ")

  # SORT
  ReAL_order <- c("Subject",

                  "correct TargetCat1 TR Compatible", "incorrect TargetCat1 TR Compatible", "correct TargetCat1 TR Incompatible", "incorrect TargetCat1 TR Incompatible",
                  "correct TargetCat2 TR Compatible", "incorrect TargetCat2 TR Compatible", "correct TargetCat2 TR Incompatible", "incorrect TargetCat2 TR Incompatible",
                  "correct AttributeCat1 TR Compatible", "incorrect AttributeCat1 TR Compatible", "correct AttributeCat1 TR Incompatible", "incorrect AttributeCat1 TR Incompatible",
                  "correct AttributeCat2 TR Compatible", "incorrect AttributeCat2 TR Compatible", "correct AttributeCat2 TR Incompatible","incorrect AttributeCat2 TR Incompatible",

                  "correct TargetCat1 TS Compatible", "incorrect TargetCat1 TS Compatible", "correct TargetCat1 TS Incompatible", "incorrect TargetCat1 TS Incompatible",
                  "correct TargetCat2 TS Compatible", "incorrect TargetCat2 TS Compatible", "correct TargetCat2 TS Incompatible", "incorrect TargetCat2 TS Incompatible",
                  "correct AttributeCat1 TS Compatible", "incorrect AttributeCat1 TS Compatible", "correct AttributeCat1 TS Incompatible", "incorrect AttributeCat1 TS Incompatible",
                  "correct AttributeCat2 TS Compatible", "incorrect AttributeCat2 TS Compatible", "correct AttributeCat2 TS Incompatible","incorrect AttributeCat2 TS Incompatible")

  ReAL_data <- ReAL_data[,..ReAL_order]



  # RECODING INDICATOR TO SUBSET
  SwitchingCosts_Comp <- tmp[is.na(TaskSwitch) == FALSE,
                             .(SwitchingCosts = mean(Accuracy[TaskSwitch == "TS" & Compatibility == "Compatible"] == 0)
                               - mean(Accuracy[TaskSwitch == "TR" & Compatibility == "Compatible"] == 0)),
                             by = "Subject"][,2]

  SwitchingCosts_Incomp <- tmp[is.na(TaskSwitch) == FALSE,
                               .(SwitchingCosts = mean(Accuracy[TaskSwitch == "TS" & Compatibility == "Incompatible"] == 0)
                                 - mean(Accuracy[TaskSwitch == "TR" & Compatibility == "Incompatible"] == 0)),
                               by = "Subject"][,2]

  ReAL_data$Re_in_comp <- SwitchingCosts_Incomp >= SwitchingCosts_Comp
  ReAL_data$Re_in_incomp <- SwitchingCosts_Incomp < SwitchingCosts_Comp
  ReAL_data$SwitchCostDif <- SwitchingCosts_Incomp - SwitchingCosts_Comp


  colnames(ReAL_data) <- c("Subject",paste0("F",sprintf("%02d",1:32)),"Re_in_comp","Re_in_incomp","SwitchCostDif")
  class(ReAL_data) <- "data.frame"

  print(
    paste(
      nrow(ReAL_data), "subjects,",
      sum(ReAL_data$Re_in_comp), "recoded in the compatible block,",
      unique(rowSums(ReAL_data[,paste0(paste0("F",sprintf("%02d",1:32)))])), "Trials."))


  return(ReAL_data)
}

#' Preparing Data for Quad Model
#'
#' This function allows you to prepare IAT data for the application of functions estimating Quad model parameters.
#' @param data Dataframe containing all relevant columns (Subject, Compatibility, StimulusCat and Correct).
#' @param Subject Name of the column containing the subject id.
#' @param Compatibility Name of the column containing the compatibility of current key mapping.
#' @param Compatible Level of Compatibility coding a compatible mapping.
#' @param Incompatible Level of Compatibility coding an incompatible mapping.
#' @param StimulusCat Name of the column containing the category of the current stimulus.
#' @param Correct Name of the column containing accuracy of responses (0 = error, 1 = correct).
#' @param TargetCat1 Level of StimulusCat coding the first target category (i.e, flower)
#' @param TargetCat2 Level of StimulusCat coding the second target category (i.e, insect)
#' @param AttributeCat1 Level of StimulusCat coding the first attribute category (i.e, pleasant)
#' @param AttributeCat2 Level of StimulusCat coding the second attribute category (i.e, unpleasant)
#' @keywords IAT mpt Quad
#' @export
#' @examples
#' prep_Quad()

prep_Quad <- function(data,
                      Subject,
                      Compatibility,
                      Compatible = 1,
                      Incompatible = 2,
                      StimulusCat,
                      Correct,
                      TargetCat1,
                      TargetCat2,
                      AttributeCat1,
                      AttributeCat2){


  if(class(data)[1] != "data.frame") stop('Data must be a data.frame!')
  if(all(c(Subject,Compatibility,StimulusCat,Correct)  %in% names(data)) == FALSE  ) stop('Check column names!')
  if(setequal(c(TargetCat1,TargetCat2,AttributeCat1,AttributeCat2), unique(data[[StimulusCat]])) == FALSE  ) stop('Check category names!')

  tmp <- data.table(Subject = data[[Subject]],
                    Compatibility = factor(data[[Compatibility]], levels = c(1,2), labels = c("Compatible","Incompatible")),
                    Accuracy = data[[Correct]])


  # tmp <- tmp[data[[Trial]]!= min(data[[Trial]])]

  tmp$StimulusCat[data[[StimulusCat]] == TargetCat1] <- "TargetCat1"
  tmp$StimulusCat[data[[StimulusCat]] == TargetCat2] <- "TargetCat2"

  tmp$StimulusCat[data[[StimulusCat]] == AttributeCat1] <- "AttributeCat1"
  tmp$StimulusCat[data[[StimulusCat]] == AttributeCat2] <- "AttributeCat2"

  Quad_data <- tmp[,.(correct = sum(Accuracy == 1),
                      incorrect = sum(Accuracy == 0)),
                   by = c("Subject", "StimulusCat", "Compatibility")]

  # -> WIDE FORMAT
  Quad_data <- dcast.data.table(Quad_data, Subject ~ StimulusCat + Compatibility, value.var = c("incorrect","correct"), sep = " ")

  # SORT
  Quad_order <- c("Subject",

                  "correct TargetCat1 Compatible", "incorrect TargetCat1 Compatible", "correct TargetCat2 Compatible", "incorrect TargetCat2 Compatible",
                  "correct AttributeCat1 Compatible", "incorrect AttributeCat1 Compatible","correct AttributeCat2 Compatible", "incorrect AttributeCat2 Compatible",
                  "correct TargetCat1 Incompatible", "incorrect TargetCat1 Incompatible", "correct TargetCat2 Incompatible", "incorrect TargetCat2 Incompatible",
                  "correct AttributeCat1 Incompatible", "incorrect AttributeCat1 Incompatible","correct AttributeCat2 Incompatible","incorrect AttributeCat2 Incompatible")

  Quad_data <- Quad_data[,..Quad_order]

  colnames(Quad_data) <- c("Subject",paste0("t",sprintf("%02d",1:16)))

  print(paste(nrow(Quad_data), "subjects,",
              unique(rowSums(Quad_data[,paste0(1:16)])), "Trials.")
  )
  class(Quad_data) <- "data.frame"


  return(Quad_data)
}
