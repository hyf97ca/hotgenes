#' Launches the shiny app for hotgenes
#'
#' A function that launches the shiny app for hotgenes.
#'
#' @return No return value, but starts a shiny server.
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @export
#' @importFrom shiny runApp

runHotgenes <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "hotgenes")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
