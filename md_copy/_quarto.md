project:
  type: website
  output-dir: docs
  render:
    - "*.qmd"
    - "!quarto_shiny.qmd"
    - "!quartoP5.qmd"
    - "!alm_hidden.qmd"
    - "!bayes*"
    - "!htw_expGPT.qmd"
    #- !htw_exam.qmd"

exclude: ["*.rds","*.tex","*.csl","*.bib","docs*","404.html","_old*","out.bib",
"*.pdf","*.R",".obsidian","Rproj","R","Bib","exp","HTW_Analysis/**","*HTW_Modelling/**","Motivated*","Manuscript*"]

# Exclude all files in HTW_Analysis folder, recursively
website:
  title: "HTW Project"
  navbar:
    pinned: false
   # background: "#EE6A24"
    left:
      - text: Analyses
        menu: 
          - discrim.qmd
          - ME_Pool.qmd
          - analysis.qmd
      - text: Simulations
        menu: 
          - simParameterRecovery.qmd
          - DeLosh97_Sim.qmd
          - SimReplications.qmd
          - alm_learning.qmd
      - text: Model Fitting
        menu: 
          - htw_exam.qmd
          - abc_htw.qmd
          - group_fits.qmd
      - text: Interactive
        menu: 
          - ojs_alm.qmd
          - quarto_shiny.qmd
          - ALM_Shiny.qmd
          - quartoP5.qmd
      - text: Misc
        menu: 
          - Task.qmd
          - model_viz.qmd
          - htw_dp.qmd
          - External.qmd
          - HTW_ToDo.qmd
          - benchmarks.qmd


format:
  html:
    theme:
      light: spacelab #[default]
      dark: cyborg
    css: ["Style/style1.css", "Style/calloutTG.css"]
    toc: true
    toc-location: left
    cold-fold: true
    cold-tool: true
    smooth-scroll: true
    #html-math-method: mathjax

execute:
  freeze: true




# website:
#   title: "HTW Project"
#   navbar:
#     pinned: false
#    # background: "#EE6A24"
#     left:
#       - text: To-Do
#         file: HTW_ToDo.qmd
#       - text: Analyses
#         file: discrim.qmd
#       - text: Simulations
#         menu: 
#           - DeLosh97_Sim.qmd
#           - SimReplications.qmd
#           - ALM_Shiny.qmd
#           - ojs_alm.qmd
#           - alm_learning.qmd
#       - text: Model Fitting
#         menu: 
#           - htw_exam.qmd
#           - benchmarks.qmd
#       - text: Task
#         file: Task.qmd
#       - text: Proposal
#         file: htw_dp.qmd
#       - External.qmd
