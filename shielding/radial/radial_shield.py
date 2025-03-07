import numpy as np
import openmc
import openmc.data
import openmc.model
import openmc.stats

power = 4.88e9
batches = 1000
# Retrieve info from statepoint
photon_energy_bins, photon_coeffs = openmc.data.dose_coefficients("photon")
neutron_energy_bins, neutron_coeffs = openmc.data.dose_coefficients("neutron")
state = openmc.StatePoint("core_statepoint.h5", autolink=False)
photon_tally = state.get_tally(id=3)
neutron_tally = state.get_tally(id=4)
# Set up model
model = openmc.model.Model()

"""Materials"""
void = openmc.Material(name="Void")
void.set_density("g/cm3", 1.0e-10)
void.add_nuclide("H3", 1)

air = openmc.Material(name="Air")
air.set_density("g/cm3", 1.1839e-3)
air.add_element("N", 0.7808)
air.add_element("O", 0.2095)
air.add_element("Ar", 0.093)
air.add_element("C", 0.0004)

water = openmc.Material(name="Water")
water.set_density("g/cm3", 1.0)
water.add_element("H", 2)
water.add_element("O", 1)

real_water = openmc.Material(name="Water")
real_water.set_density("g/cm3", 1.0)
real_water.add_element("H", 2)
real_water.add_element("O", 1)
real_water.add_s_alpha_beta("c_H_in_H2O")

concrete = openmc.Material(name="Concrete")
concrete.set_density("g/cm3", 2.95)
concrete.add_element("Si", 0.1699779857, "wo")
concrete.add_element("Al", 0.008564609817, "wo")
concrete.add_element("Fe", 0.2443627177, "wo")
concrete.add_element("Ca", 0.04064262112, "wo")
concrete.add_element("Mg", 0.07204427895, "wo")
concrete.add_element("S", 0.001433048006, "wo")
concrete.add_element("Cl", 0.0001417671513, "wo")
concrete.add_element("Na", 0.001946852724, "wo")
concrete.add_element("K", 0.002117051692, "wo")
concrete.add_element("Ti", 0.0003337333589, "wo")
concrete.add_element("P", 0.002345614791, "wo")
concrete.add_element("Mn", 0.001388663595, "wo")
concrete.add_element("O", 0.4468310563, "wo")
concrete.add_element("H", 0.00767, "wo")

steel = openmc.Material(name="Steel")
steel.set_density("g/cm3", 8.0)
steel.add_element("C", 0.07, "wo")
steel.add_element("Cr", 18.50, "wo")
steel.add_element("Mn", 0.2, "wo")
steel.add_element("Ni", 9.25, "wo")
steel.add_element("P", 0.045, "wo")
steel.add_element("S", 0.03, "wo")
steel.add_element("Si", 1.0, "wo")
steel.add_element("Fe", 70.905, "wo")

reactor_compartment = openmc.Material.mix_materials(
    [water, air, steel], (0.35, 0.64, 0.01), "vo", name="Reactor Compartment"
)
reactor_compartment.add_s_alpha_beta("c_H_in_H2O")

primary_shield = openmc.Material(name="Primary Shield")
primary_shield.set_density("g/cm3", 11.348)
primary_shield.add_element("Pb", 1)

secondary_shield = openmc.Material.mix_materials(
    [concrete], [1], name="Secondary Shield"
)
secondary_shield.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials(
    [
        void,
        air,
        real_water,
        steel,
        reactor_compartment,
        primary_shield,
        secondary_shield,
    ]
)
model.materials = materials

"""Geometry"""
reactor_vessel_or = openmc.ZCylinder(x0=0, y0=0, r=249.0, name="Reactor Outer Radius")
inner_lead_or = openmc.ZCylinder(x0=0, y0=0, r=255.0, name="Lead Lining")
pswt_or = openmc.ZCylinder(x0=0, y0=0, r=310.0, name="PSWT Outer Radius")
primary_shield_or = openmc.ZCylinder(
    x0=0, y0=0, r=332.0, name="Primary Shield Outer Radius"
)
steel_lining_ir = openmc.ZCylinder(x0=0, y0=0, r=1867)
secondary_shield_ir = openmc.ZCylinder(x0=0, y0=0, r=1870)
secondary_shield_or = openmc.ZCylinder(x0=0, y0=0, r=1966)
max_radius = openmc.ZCylinder(x0=0, y0=0, r=2000, boundary_type="vacuum")

reactor_cell = openmc.Cell(fill=void, region=-reactor_vessel_or)
lead_lining = openmc.Cell(
    fill=primary_shield, region=+reactor_vessel_or & -inner_lead_or
)
pswt_cell = openmc.Cell(fill=real_water, region=+inner_lead_or & -pswt_or)
primary_cell = openmc.Cell(fill=primary_shield, region=+pswt_or & -primary_shield_or)
rc_cell = openmc.Cell(
    fill=air,
    region=+primary_shield_or & -steel_lining_ir,
)
steel_lining_cell = openmc.Cell(
    fill=steel, region=+steel_lining_ir & -secondary_shield_ir
)
secondary_cell = openmc.Cell(
    fill=secondary_shield, region=+secondary_shield_ir & -secondary_shield_or
)
person_cell = openmc.Cell(fill=real_water, region=+secondary_shield_or & -max_radius)

root = openmc.Universe(name="root universe")
root.add_cells(
    (
        reactor_cell,
        lead_lining,
        pswt_cell,
        primary_cell,
        rc_cell,
        steel_lining_cell,
        secondary_cell,
        person_cell,
    )
)

geometry = openmc.Geometry(root)
model.geometry = geometry

"""Settings"""
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.batches = batches
settings.particles = 100000
# Photon source
photon_source = openmc.IndependentSource()
photon_source.space = openmc.stats.Point((0, 0, 0))
photon_source.strength = 0.726065
photon_source.particle = "photon"
photon_source.energy = openmc.EnergyFilter(photon_energy_bins).get_tabular(
    photon_tally.mean[:, 0, 0]
)
# Neutron source
neutron_source = openmc.IndependentSource()
neutron_source.space = openmc.stats.Point((0, 0, 0))
neutron_source.strength = 0.273935
neutron_source.energy = openmc.EnergyFilter(neutron_energy_bins).get_tabular(
    neutron_tally.mean[:, 0, 0]
)
nitrogen_source = openmc.IndependentSource()
nitrogen_source.space = openmc.stats.Point((0, 0, 0))
nitrogen_source.strength = 5.98615e-3
nitrogen_source.particle = "photon"
nitrogen_source.energy = openmc.stats.Discrete([6.128e6, 7.115e6], [67.9, 4.9])
settings.source = [photon_source, neutron_source, nitrogen_source]
settings.photon_transport = True
model.settings = settings

"""Plots"""
plot1 = openmc.Plot()
plot1.origin = (0, 0, 0)
plot1.width = (4000, 4000)
plot1.pixels = (4000, 4000)
plot1.color_by = "material"
model.plots = [plot1]

"""Tallies"""
out_cell_filter = openmc.CellFilter(person_cell)
photon_energy_filter = openmc.EnergyFunctionFilter(photon_energy_bins, photon_coeffs)
photon_particle_filter = openmc.ParticleFilter("photon")
neutron_energy_filter = openmc.EnergyFunctionFilter(neutron_energy_bins, neutron_coeffs)
neutron_particle_filter = openmc.ParticleFilter("neutron")
tally1 = openmc.Tally(1)
tally1.filters = [out_cell_filter, photon_energy_filter, photon_particle_filter]
tally1.scores = ["flux"]
tally2 = openmc.Tally(2)
tally2.filters = [out_cell_filter, neutron_energy_filter, neutron_particle_filter]
tally2.scores = ["flux"]
model.tallies = [tally1, tally2]

"""Execution"""

model.plot_geometry()
model.run(export_model_xml=False)
results = openmc.StatePoint("statepoint.{}.h5".format(batches))
photon_counts = results.get_tally(id=1)
neutron_counts = results.get_tally(id=2)
heat_tally = state.get_tally(id=2)
heating = heat_tally.mean[0, 0, 0] * 1.602e-19
source_strength = power / 4 / heating
circum = np.pi * (1965.96 + 2000)
volume = np.pi * (2000**2 - 1965.96**2)
out_photon_dose = photon_counts.mean[:, 0, 0] * source_strength / volume / circum * 0.5
out_neutron_dose = (
    neutron_counts.mean[:, 0, 0] * source_strength / volume / circum * 0.5
)
out_total_dose = np.sum(out_photon_dose) + np.sum(out_neutron_dose)
print("Outside photon dose rate: %g uSv/h" % (np.sum(out_photon_dose) * 3600 * 1e-6))
print("Outside neutron dose rate: %g uSv/h" % (np.sum(out_neutron_dose) * 3600 * 1e-6))
print("Outside total dose rate: %g uSv/h" % (out_total_dose * 3600 * 1e-6))
