## Code to make PCR reactions using a pre-mixed master-mix, distribute it between wells
## of a 96-well plate, add primers and add a given number of template samples.
## each reaction will be performed in replicate in order to perform a gradient 
## of annealing temps

##Protocol##
#1 - make master mix in trough (contains premixed MM, primers, water and DNA) using p200
#2 - transfer master mix to wells of a 96 wells plate using p300 multi-channel

###INPUT### PCR variables
num_replicates = 5
num_templates = 10

total_pcr_volume = 50
master_mix_volume = 20
template_volume = 2
primer_volume = 2.5

from opentrons import robot, containers, instruments

#Define containers - source_tube rack = cold block
mix_trough = containers.load('96-deep-well', 'B2', 'mix-trough')
source_trough = containers.load('trough-12row', 'B1', 'source-trough')
trash = containers.load('trash-box', 'A3')

p200rack = containers.load('tiprack-200ul', 'C3', 'p200-rack')
mc_p300rack = containers.load('tiprack-200ul', 'E1', 'mc-p300_rack')
mc_p300rack2 = containers.load('tiprack-200ul', 'E2', 'mc-p300_rack2')

#Create 8x12 PCR plate
containers.create(
	'PCR-plate',
	grid=(8,12),
	spacing=(9,9),
	diameter=5,
	depth=19.5
)

output = containers.load('PCR-plate', 'D1', 'output')

#Create 3x6 2ml tube rack for DNA samples
containers.create(
	'3x6-tube-rack-2ml',
	grid=(3,6),
	spacing=(19.5,19.5),
	diameter=9.5,
	depth=40
)

template_rack = containers.load('3x6-tube-rack-2ml', 'B3', 'template-rack')

#Define single channel pipette
p200 = instruments.Pipette(
    trash_container=trash,
    tip_racks=[p200rack],
    min_volume=20,
    max_volume=200,
    axis="b"
)

#Define multichannel pipette
p300_mc = instruments.Pipette(
	trash_container=trash,
	tip_racks=[mc_p300rack, mc_p300rack2],
	min_volume=50,
	max_volume=300,
	axis="a",
	channels=8
)

#Define DNA volumes
import math
num_templates_nearest8 = int(math.ceil(num_templates / 8)) * 8
template_volumes = [template_volume] * num_templates
num_pcr_samples = len(template_volumes)

num_templates2 = num_templates_nearest8 / 8
num_templates3 = int(round(num_templates2))
num_rows = [template_volume] * num_templates3
num_rows2 = len(num_rows)

#Define locations of PCR components
water_source = source_trough.wells('A1')
pcr_mix_source = source_trough.wells('A2')
F_primer_source = source_trough.wells('A3')
R_primer_source = source_trough.wells('A4')
template_sources = template_rack.wells('A1', length=num_pcr_samples)
complete_mix_sources = mix_trough.rows('1', length=num_pcr_samples)

#Define components of master mix
total_mm_volume = (total_pcr_volume * (num_replicates+1))
mm_pcr_mix_volume = (master_mix_volume * (num_replicates+1))
mm_water_volume = ((total_pcr_volume - master_mix_volume - (2*primer_volume) - (template_volume)) * (num_replicates+1))
mm_primer_volume = (primer_volume * (num_replicates+1))
mm_primer_volume = (primer_volume * (num_replicates+1))
mm_template_volume = (template_volume * (num_replicates+1))

#Define PCR-plate wells
PCR_plate_output_wells = [well.bottom(1) for well in output.rows('1', length=num_replicates)]

#Define mix wells
mix_wells1 = [well.bottom() for well in mix_trough.rows('1', length=num_rows2)]
mix_wells2 = [well.bottom() for well in mix_trough.wells('A1', length=num_pcr_samples)]

#Make master mix
p300_mc.distribute(
	mm_water_volume, water_source, mix_wells1, blow_out=True)

p300_mc.transfer(
	mm_pcr_mix_volume, pcr_mix_source, mix_wells1, mix_after=(2, 50), blow_out=True,
	new_tip='always')

p300_mc.transfer(
	mm_primer_volume, F_primer_source, mix_wells1, mix_after=(2, 50), blow_out=True,
	new_tip='always')

p300_mc.transfer(
	mm_primer_volume, R_primer_source, mix_wells1, mix_after=(2, 50), blow_out=True,
	new_tip='always')

p200.transfer(
	mm_template_volume, template_sources, mix_wells2, mix_after=(2, 50), blow_out=True, 
	new_tip='always')

p300_mc.pick_up_tip()
p300_mc.mix(
	5,300, mix_trough.wells('A1', length=num_pcr_samples))

#Transfer master mixes to 96-well plate
p300_mc.distribute(
	total_pcr_volume, complete_mix_sources, PCR_plate_output_wells, 
	disposal_vol=10, dispense_speed=300)
