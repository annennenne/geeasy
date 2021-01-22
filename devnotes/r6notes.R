library(R6)

Person <- R6Class("Person",
                  public = list(
                    name = NA,
                    hair = NA,
                    initialize = function(name, hair) {
                      if (!missing(name)) self$name <- name
                      if (!missing(hair)) self$hair <- hair
                      self$greet()
                    },
                    set_hair = function(val) {
                      self$hair <- val
                    },
                    greet = function() {
                      cat(paste0("Hello, my name is ", self$name, ".\n"))
                    }
                  )
)


joe <- Person$new(name = "joe")
ann <- Person$new(name = "ann")

joe2 <- joe$clone()
joe2$set_hair("blue")


joe
kylie <- Person$new("kylie")
kylie$greet()
