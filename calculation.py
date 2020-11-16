# Imports
from astropy import units as u
import numpy as np


# Constants
NA = 6.02214076e+23 / u.mol
TOTAL_MASS = 6.25e21 * u.gram  # Coustenis
DEMAND_YEAST = 3.2 * 1e3*u.joule / u.gram
GENTIME_MAX = 1 * u.year


# Classes
class Molecule:
    def __init__(self, name, mole_frac, molar_mass):
        '''
        mole_frac  : mole fraction of the molecule in the mixture (decimal, not %)
        molar_mass : molar mass of the molecule in g/mol
        '''
        self.name = name
        self.mole_frac = mole_frac
        self.molar_mass = molar_mass * u.gram/u.mol

class Reaction:
    def __init__(self, name, methane_produced, energy_per_mol):
        '''
        For reactions that produce methane only.
        methane_produced : number of methane produced as in the chemical formula
        energy_per_mol   : the energy released in kJ/mol
        '''
        self.name = name
        self.methane_produced = methane_produced
        self.energy_per_mol = energy_per_mol * 1e3 * u.joule / u.mol


# Globals used in functions
N2 = Molecule('N2', 0.95, 28.0134)
CH4 = Molecule('CH4', 0.049, 16.04246)
H2 = Molecule('H2', 0.001, 2.01588)
ATM_MOLECULES = [N2, CH4, H2]
R1 = Reaction('1', 2, 334)
R2 = Reaction('2', 2, 57)
ALL_REACTIONS = [R1, R2]


# Functions
def mass_in_mixture(molecule):
    '''
    Total mass of the molecule in the mixture
    molecule : Molecule object
    '''
    den = 0
    for m in ATM_MOLECULES:
        den += m.mole_frac * m.molar_mass
    num = molecule.mole_frac * molecule.molar_mass
    return num / den * TOTAL_MASS

def moles(molecule):
    '''
    Number of moles of the molecule (Molecule object) in the mixture
    '''
    return mass_in_mixture(molecule) / molecule.molar_mass

def production_rate_mol(molecule, tau):
    '''
    Production rate to balance out exponential decay, in mol/yr
    '''
    N0 = moles(molecule)
    t = 1 * u.year
    return N0 / tau * np.exp(- t / tau)

def get_reaction_ratios(areas, perimeters, peri_width):
    '''
    Given the inner areas and perimeter lengths of multiple lakes,
    gets the fractions of perimeter area and inner area over the total area,
    also returns the total area.
    areas      : array of the (inner) surface areas of lakes
    perimeters : array of the perimeter lengths of lakes
    peri_width : ribbon width of the perimeters
    '''
    perimeters = perimeters * peri_width
    surface_total = areas.sum() + perimeters.sum()
    acetylene = perimeters.sum() / surface_total
    ethane = areas.sum() / surface_total
    return acetylene, ethane, surface_total

def production_rate_kJ(reactions, reaction_ratios, tau):
    '''
    Converts production rate from mol/yr to kJ/yr
    reactions       : list of Reaction objects
    reaction_ratios : list of ratios for each reaction (percent taking place in 
                      unit area w.r.t. all reactions)
    '''
    mol_per_yr = production_rate_mol(CH4, tau)
    mol_test = 0
    for i in range(len(reactions)):
        mol_test += reaction_ratios[i] * reactions[i].methane_produced
    energy = 0
    for i, reaction in enumerate(reactions):
        frac = reaction_ratios[i] * reaction.methane_produced / mol_test
        energy += mol_per_yr * frac / reaction.methane_produced * reaction.energy_per_mol
    return energy

def get_demand(energy_per_mass, generation_time):
    '''
    Gets the microbe's energy demand in (J/g/yr)
    energy_per_mass : energy it takes to build a unit mass of microbe
    generation_time : generation time of the microbe species
    '''
    return energy_per_mass / generation_time

def biomass_density(energy, demand, area):
    '''
    Gets the biomass surface density (g/m2)
    energy : energy produced per year (J/yr)
    demand : energy demand of the microbes (J/g/yr)
    area   : populated area (m2)
    '''
    return energy / demand / area

def get_dry_mass(choice):
    '''
    Returns dry mass of microbe based on acrylonitrile vesicles or Earth average
    '''
    if choice == 'acryl':
        # Palmer et al 2017
        acryl_ligeia = 1e14 * u.kilogram
        density_ligeia = (3e7 / u.cm**3).to(1/u.m**3)
        volume_ligeia = 14000 * u.kilometer**3

        dry_mass_acryl = acryl_ligeia / (density_ligeia * volume_ligeia)
        return dry_mass_acryl.to(u.gram)

    if choice == 'earth':
        return 2e-14 * u.gram

def microbe_density(energy, demand, area, dry_mass):
    '''
    Calculates microbes number density
    energy   : energy produced per year (J/yr)
    demand   : energy demand of the microbes (J/g/yr)
    area     : populated area (m2)
    dry_mass : mean dry mass of microbes (g)
    '''
    return biomass_density(energy, demand, area) / dry_mass


# Main routine
def main(tau, peri_width, energy_per_mass, generation_time, choice):

    print(f"Moles of methane in atmosphere = {moles(CH4):.2e}")
    print(f"tau = {tau}")
    print(f"Moles of methane produced per year = {production_rate_mol(CH4, tau):.2e}")

    lake_data = np.loadtxt('shoreline.txt', unpack = True)
    areas = lake_data[3] * u.kilometer**2
    perimeters = lake_data[4] * u.kilometer
    R1_frac, R2_frac, POP_SURFACE = get_reaction_ratios(areas, perimeters, peri_width)
    print(f"Reaction 1 fraction = {R1_frac:.3f}")
    print(f"Reaction 2 fraction = {R2_frac:.3f}")
    print(f"Total populated surface = {POP_SURFACE}")

    energy = production_rate_kJ([R1, R2], [R1_frac, R2_frac], tau)
    energy_max = production_rate_kJ([R1], [1], tau)
    energy_min = production_rate_kJ([R2], [1], tau)
    print(f"Minimum energy (100% Reaction 2) = {energy_min:.2e}")
    print(f"Maximum energy (100% Reaction 1) = {energy_max:.2e}")
    print(f"According to fractions, energy = {energy:.2e}")

    demand = get_demand(energy_per_mass, generation_time)
    grams_per_area = biomass_density(energy, demand, POP_SURFACE).to(u.gram/u.meter**2)
    print(f"biomass density = {grams_per_area:.2e}")

    dry_mass = get_dry_mass(choice)
    print(f"with chosen dry mass as {choice} i.e. {dry_mass:.2e}")
    microbe_per_area = microbe_density(energy, demand, POP_SURFACE, dry_mass)
    print(f"number density = {microbe_per_area.to(1/u.meter**2):.2e}")
    print()
    

# Run main with different params
tau_kras = 2.13e7 * u.year
tau_waite = ((moles(CH4) * NA) / (5e27 / u.second)).to(u.year)
tau_miller = (moles(CH4) / (7e9 * u.kilogram / u.year / CH4.molar_mass)).to(u.year)

for i in [tau_kras, tau_waite, tau_miller]:
    for j in ['acryl', 'earth']:
        main(i, 0.5 * u.kilometer, DEMAND_YEAST / 2, GENTIME_MAX / 2, j)
