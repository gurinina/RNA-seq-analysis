# Render the book in GitBook format
bookdown::render_book("index.Rmd", "bookdown::gitbook")

# Open the resulting GitBook in the Viewer pane
rstudioapi::viewer("_book/index.html")
