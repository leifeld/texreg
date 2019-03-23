#' Publish praise about \pkg{texreg}
#'
#' Publish praise about \pkg{texreg} to help the developers demonstrate impact.
#'
#' You can use this function to praise the \pkg{texreg} package. Funders and
#' academic employers are increasingly interested in seeing evidence for the
#' impact academic research generates. For software, such as \pkg{texreg}, this
#' is very hard to accomplish because the developers are usually disconnected
#' from the users. The consequence is that incentives for developing
#' packages like these are diminishing the more the funders and employers
#' require evidence of impact on society, firms, or policy makers.
#'
#' The \code{\link{praise}} function is our attempt at rectifying the situation.
#' With this function, you can provide positive feedback to the developers. The
#' praise is saved to a database on the web server of the package maintainer and
#' subsequently displayed at \url{http://www.philipleifeld.com/texreg-praise}
#' for other users, funders, and employers to view. This will also enable the
#' package authors to compile reports about how \pkg{texreg} is used by academic
#' and non-academic users to increase their productivity and work quality, for
#' example in the form of an impact case study for the 2021 UK Research
#' Excellence Framework (REF).
#'
#' We need many positive examples of how \pkg{texreg} has an impact on your
#' work. We are especially interested in non-academic users, but welcome
#' feedback from anyone. So please contribute by using the praise function! Tell
#' us how cool this package is and how it has changed your work!
#'
#' The minimal information we require from you is whether you are an academic or
#' non-academic user, the name of your organization, and some free-form praise
#' (of a general nature, or about how it makes you more productive, or about how
#' it increases the quality of your work or reporting). But there are some
#' additional fields. While we are happy with the basic information, of course
#' we will be happier if we also know your name, how to contact you, what kinds
#' of models you work with, and some other details. Your choice!
#'
#' Please note that by using the \code{\link{praise}} function you agree that
#' the information you provide through the function, including your location, is
#' stored online in a database, displayed on the website of the package author,
#' and used in reports to funders, employers etc. (This is the whole purpose of
#' it.) You can contact the package maintainer any time to have your praise
#' removed.
#'
#' @param academic_user Should be \code{TRUE} if you are at a university or
#'   public research institute. Should be \code{FALSE} if you are a private
#'   user, for example you are using \pkg{texreg} in your work for a firm, NGO,
#'   association, government department, as an individual user etc. We
#'   particularly need praise from non-academic users to demonstrate societal
#'   impact, but we can also make the case for academic usage to generate impact
#'   indirectly.
#' @param organization Please tell us the name of the organization for which you
#'   are using \pkg{texreg}. If we can show that the package is being employed
#'   in a number of different settings, this will help us demonstrate impact.
#' @param general_praise Use this argument to provide general praise, for
#'   example about the way it was designed, the user support you have received,
#'   or just how much you enjoy using it. While this is useful, however, we
#'   would be even more interested in receiving statements in how \pkg{texreg}
#'   makes you more productive (in the \link{increase_productivity} argument) or
#'   how it increases the quality of your work or your reports (through the
#'   \link{increase_quality} function). Note: you need to provide at least one
#'   of these three free-form text arguments.
#' @param increase_productivity This is one of the fields we are most interested
#'   in. Please use this field to tell us how \pkg{texreg} is making you more
#'   productive. For example, does it speed up writing your articles or research
#'   reports? Does it enable you to skip manual work like copy and paste of your
#'   results into your reports, or to avoid fiddling with table formatting? How
#'   much time has it saved you so far? Are there any other benefits in terms of
#'   productivity you can think of? Note: you need to provide feedback using at
#'   least one of the three free-form arguments (\code{general_praise},
#'   \code{increase_productivity}, or \code{increase_quality}).
#' @param increase_quality This is one of the fields we are most interested in.
#'   Please use this argument to tell us how \pkg{texreg} increases the quality
#'   of your work or the quality of your reporting. For example, does the
#'   package generate tables that look more professional than the tables you
#'   used to create manually? Are you using \link{screenreg} to improve your
#'   workflow by understanding better how the results of multiple models
#'   compare? Are you using \link{plotreg} to visualize and present your
#'   statistical results in a more effective way? Can you think of any other
#'   ways in which \pkg{texreg} is helping you? Note: you need to provide
#'   feedback using at least one of the three free-form arguments
#'   (\code{general_praise}, \code{increase_productivity}, or
#'   \code{increase_quality}).
#' @param start_using (Optional) When did you start using \pkg{texreg}? We are
#'   interested in the approximate time or year as a free-form text argument,
#'   for example \code{"back in 2013 when the JSS article came out"}.
#' @param where_learn (Optional) Where or how did you learn about the
#'   \pkg{texreg} package?
#' @param name (Optional) We would be delighted to to know who you are. After
#'   all, we can quote you much more effectively if we can tell the funders and
#'   employers who provided this praise! If possible, include your title.
#' @param contact_details (Optional) Tell us how we can contact you in case we
#'   would benefit from additional information. This might help us further down
#'   the road in compiling an impact case study or a similar report. Don't
#'   worry, this information will not be displayed on the website!
#' @param models (Optional) Which kinds of statistical models do you use in your
#'   work? For example, \code{"Mostly linear models, but also lme4 and ergm."}.
#' @param num_users (Optional) How many other \pkg{texreg} users do you know? In
#'   particular, if you are a non-academic user, would you mind telling us how
#'   many other non-academic users you are aware of and how many of them are in
#'   your organization? The more we know, the more convincing our evidence base
#'   will be. This argument accepts \code{numeric} values or more detailed
#'   responses as a \code{character} object.
#' @return If everything works well, no output is returned. If the submission of
#'   the praise to the maintainer fails, a \code{response} object (as defined in
#'   the \pkg{httr} package) will be returned. Should you have any problems, do
#'   feel free to e-mail your praise to the package maintainer directly.
#'
#' @author Philip Leifeld
#' @aliases praise, feedback
#'
#' @examples
#' \dontrun{
#' praise(academic_user = TRUE,
#'        organization = "University of Happy Tables",
#'        increase_quality = "Man I've never had such pretty tables!")
#' }
#'
#' @import httr
#' @export
praise <- function(academic_user,
                   organization,
                   general_praise = NULL,
                   increase_productivity = NULL,
                   increase_quality = NULL,
                   start_using = NULL,
                   where_learn = NULL,
                   name = NULL,
                   contact_details = NULL,
                   models = NULL,
                   num_users = NULL) {

  # process 'academic_user' argument
  if (is.null(academic_user) || is.na(academic_user) || !is.logical(academic_user) || length(academic_user) != 1) {
    stop("'academic_user' must be TRUE if you are at a university or research institute and FALSE otherwise.")
  }
  user_type <- "N/A"
  if (isTRUE(academic_user)) {
    user_type <- "academic"
  } else {
    user_type <- "non-academic"
  }

  # process 'organization' argument
  if (is.null(organization) || is.na(organization) || !is.character(organization) || length(organization) != 1) {
    stop("'organization' must be the name of your organization. You may write 'private' if your texreg use is restricted to private purposes.")
  }

  # process 'general praise', 'increase_productivity', and 'increase_quality' arguments
  if (is.null(general_praise) && is.null(increase_productivity) && is.null(increase_quality)) {
    stop(paste0("At least one of the arguments 'general_praise', 'increase_productivity', or 'increase_quality' must be provided.\n",
                "- 'general_praise': Any praise you want to provide.\n",
                "- 'increase_productivity': Tell us how texreg makes you more productive.\n",
                "- 'increase_quality': Tell us how texreg increases the quality of your work or reporting."))
  }
  if (!is.null(general_praise) && (is.na(general_praise) || !is.character(general_praise) || length(general_praise) != 1)) {
    stop("'general_praise' is a free-form text argument where you can provide any kind of praise.")
  }
  if (is.null(general_praise)) {
    general_praise <- ""
  }
  if (!is.null(increase_productivity) && (is.na(increase_productivity) || !is.character(increase_productivity) || length(increase_productivity) != 1)) {
    stop("'increase_productivity' argument: Tell us how texreg makes you more productive.")
  }
  if (is.null(increase_productivity)) {
    increase_productivity <- ""
  }
  if (!is.null(increase_quality) && (is.na(increase_quality) || !is.character(increase_quality) || length(increase_quality) != 1)) {
    stop("'increase_quality' argument: Tell us how texreg increases the quality of your work or reporting.")
  }
  if (is.null(increase_quality)) {
    increase_quality <- ""
  }

  # process 'start_using' argument
  if (!is.null(start_using) && (is.na(start_using) || !is.character(start_using) || length(start_using) != 1)) {
    stop("Optional 'start_using' argument: Tell us when you started using texreg (as a character string).")
  }
  if (is.null(start_using)) {
    start_using <- ""
  }

  # process 'where_learn' argument
  if (!is.null(where_learn) && (is.na(where_learn) || !is.character(where_learn) || length(where_learn) != 1)) {
    stop("Optional 'where_learn' argument: Tell us where or how you learned about texreg (as a character string).")
  }
  if (is.null(where_learn)) {
    where_learn <- ""
  }

  # process 'name' argument
  if (!is.null(name) && (is.na(name) || !is.character(name) || length(name) != 1)) {
    stop("Optional 'name' argument: You don't have to leave your name, but we'd be happy if you did!")
  }
  if (is.null(contact_details)) {
    contact_details <- ""
  }

  # process 'contact_details' argument
  if (!is.null(contact_details) && (is.na(contact_details) || !is.character(contact_details) || length(contact_details) != 1)) {
    stop("Optional 'contact_details' argument: How can we reach you if we need to ask for further testimony for our reporting purposes?")
  }
  if (is.null(contact_details)) {
    contact_details <- ""
  }

  # process 'models' argument
  if (!is.null(models) && (is.na(models) || !is.character(models) || length(models) != 1)) {
    stop("Optional 'models' argument: What kinds of statistical models do you use in your work most often?")
  }
  if (is.null(models)) {
    models <- ""
  }

  # process 'num_users' argument
  if (!is.null(num_users) && (is.na(num_users) || length(num_users) != 1 || (!is.numeric(num_users) && !is.character(num_users)))) {
    stop("Optional 'num_users' argument: How many other texreg users do you know approximately? Supply a number or some text.")
  }
  if (is.null(num_users)) {
    num_users <- ""
  } else {
    num_users <- as.character(num_users)
  }

  # collate all information and submit
  query <- c("user_type" = user_type,
             "organization" = organization,
             "general_praise" = general_praise,
             "increase_productivity" = increase_productivity,
             "increase_quality" = increase_quality,
             "start_using" = start_using,
             "where_learn" = where_learn,
             "name" = name,
             "contact_details" = contact_details,
             "models" = models,
             "num_users" = num_users)
  url <- "http://www.webassistent.net/?"
  tryCatch({
      for (i in 1:length(query)) {
        if (i > 1) {
          url <- paste0(url, "&")
        }
        url <- paste0(url, names(query)[i], "=", URLencode(query[i], reserved = TRUE))
      }
    },
    error = function(e) message("URL could not be encoded. Did you use any exotic characters?")
  )
  tryCatch(response <- POST(url = url),
           error = function(e) message("Praise could not be sent. Please check your internet connection."))

  # user feedback
  if (grepl("<p>User type: (non-academic)|(academic)</p>", as.character(httr::content(response)))) {
    message("Success! Thank you so much for providing feedback.")
    message("This will be valuable in demonstrating the impact of texreg to funders.")
    message("You can view all feedback at http://www.philipleifeld.com/texreg-praise")
  } else {
    message("Oops, something went wrong...")
    message("Returning website response now as an httr 'response' object.")
    return(response)
  }
}