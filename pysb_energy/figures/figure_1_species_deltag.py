from ..model.example1 import model
import pysb.pattern
import pysb.bng

pysb.bng.generate_equations(model)

# TODO Calculate number of matches for each energypattern in each species and
# sum up deltaG.
