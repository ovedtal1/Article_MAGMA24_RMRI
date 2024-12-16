
# We assume TE=10 ms, slice thickness of 50 mm. 90° phase shift of first rf_ex with ADC phase folowing. Readout spoilers and linear reordering.

# Fig 3 correct TSE, 90°-180°

base_resolution= 64 # @param {type: "slider", min: 2, max: 112,step:2}
TE_ms=10 # @param {type: "slider", min: 0.0, max: 200.0}
TE=TE_ms*1e-3
TI_s=0 # @param {type: "slider", min: 0.0, max: 10.0, step: 0.1}
Excitation_FA=90 # @param {type: "slider", min: 10, max: 270}
Excitation_phase=90 # @param {type: "slider", min: 0, max: 270,step:10}
ADC_phase='same as rfex' # @param ['same as rfex','alternating']
Refocusing_FA=180 # @param {type: "slider", min: 10, max: 270}
Refocusing_phase=0 # @param {type: "slider", min: 0, max: 270,step:10}
r_spoil =1 # @param {type: "slider", min: 0, max: 3}
PEtype = 'linear' # @param ['centric', 'linear']
PE_grad_on=True # @param {type: "boolean"}
RO_grad_on=True # @param {type: "boolean"}

import numpy as np
import matplotlib.pyplot as plt
import torch
import pypulseq as pp
import MRzeroCore as mr0

# %% S1. SETUP sys
system = pp.Opts(
    max_grad=28, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
    rf_ringdown_time=20e-6, rf_dead_time=100e-6,
    adc_dead_time=20e-6, grad_raster_time=10e-6)

seq = pp.Sequence(system)
# Define FOV and resolution
fov = 200e-3
slice_thickness = 50e-3

Nread  = base_resolution  # frequency encoding steps/samples
Nphase = base_resolution  # phase encoding steps/samples


# Define rf events
rf1, gz1, gzr1 = pp.make_sinc_pulse(
    flip_angle=Excitation_FA * np.pi / 180, phase_offset=Excitation_phase * np.pi / 180, duration=1e-3,
    slice_thickness=slice_thickness, apodization=0.5, time_bw_product=4,
    system=system, return_gz=True)

rf2, gz2, _ = pp.make_sinc_pulse(
    flip_angle=Refocusing_FA* np.pi / 180, phase_offset=Refocusing_phase * np.pi / 180, duration=1e-3,
    slice_thickness=slice_thickness, apodization=0.5, time_bw_product=4,
    system=system, return_gz=True)

dwell=50e-6*2

G_flag=(int(RO_grad_on),int(PE_grad_on))  # gradient flag (read,PE), if (0,0) all gradients are 0, for (1,0) PE is off

# Define other gradients and ADC events
gx = pp.make_trapezoid(channel='x', rise_time = 0.5*dwell, flat_area=Nread / fov*G_flag[0], flat_time=Nread*dwell, system=system)
adc = pp.make_adc(num_samples=Nread, duration=Nread*dwell, phase_offset=rf1.phase_offset, delay=0*gx.rise_time, system=system)
gx_pre0 = pp.make_trapezoid(channel='x', area=+((1.0 + r_spoil) * gx.area / 2) , duration=1.5e-3, system=system)
gx_prewinder = pp.make_trapezoid(channel='x', area=+(r_spoil * gx.area / 2), duration=1e-3, system=system)
gp = pp.make_trapezoid(channel='y', area=0 / fov, duration=1e-3, system=system)
rf_prep = pp.make_block_pulse(flip_angle=180 * np.pi / 180, duration=1e-3, system=system)

if PEtype == 'centric':
      phenc = np.asarray([i // 2 if i % 2 == 0 else -(i + 1) // 2 for i in range(Nphase)]) / fov
else:
      phenc = np.arange(-Nphase // 2, Nphase // 2) / fov

# the minimal TE is given by one full period form ref pulse to ref pulse, thus gz2+gx+2*gp
minTE2=(pp.calc_duration(gz2) +pp.calc_duration(gx) + 2*pp.calc_duration(gp))/2

minTE2=round(minTE2/10e-5)*10e-5


# to realize longer TE,  we introduce a TEdelay that is added before and afetr the encoding period
TEd=round(max(0, (TE/2-minTE2))/10e-5)*10e-5  # round to raster time

if TEd==0:
  print('echo time set to minTE [ms]', 2*(minTE2 +TEd)*1000)
else:
  print(' TE [ms]', 2*(minTE2 +TEd)*1000)


# FLAIR
if TI_s>0:
  seq.add_block(rf_prep)
  seq.add_block(gx_pre0)
  seq.add_block(pp.make_delay(TI_s))


seq.add_block(rf1,gz1)
seq.add_block(gx_pre0,gzr1)

# last timing step is to add TE/2 also between excitation and first ref pulse
# from pulse top to pulse top we have already played out one full rf and gx_pre0, thus we substract these from TE/2
seq.add_block(pp.make_delay((minTE2 +TEd ) - pp.calc_duration(gz1)-pp.calc_duration(gx_pre0)))

for ii, encoding in enumerate(phenc):  # e.g. -64:63
    gp  = pp.make_trapezoid(channel='y', area=+encoding*G_flag[1], duration=1e-3, system=system)
    gp_ = pp.make_trapezoid(channel='y', area=-encoding*G_flag[1], duration=1e-3, system=system)

    seq.add_block(rf2,gz2)
    seq.add_block(pp.make_delay(TEd)) # TE delay
    seq.add_block(gx_prewinder, gp)
    if ADC_phase=='alternating':
      adc.phase_offset+=np.pi
    seq.add_block(adc, gx)
    seq.add_block(gx_prewinder, gp_)
    seq.add_block(pp.make_delay(TEd)) # TE delay

# %% S2. CHECK, PLOT and WRITE the sequence  as .seq
# Check whether the timing of the sequence is correct
ok, error_report = seq.check_timing()
if ok:
    print('Timing check passed successfully')
else:
    print('Timing check failed. Error listing follows:')
    [print(e) for e in error_report]


# %% S3 plot

seq.plot(plot_now=False,time_range=(TI_s,TI_s+2.5*TE))

# Get figure handles
fig_handles = plt.get_fignums()


# Iterate and save each figure
for fig_num in fig_handles:
    plt.figure(fig_num)
    plt.savefig(f'TSE_{fig_num}.png', format='png')
plt.close()

import subprocess
subprocess.run(["inkscape", "--export-filename=Fig_3_TSE_Pulseq.pdf", "Fig_3_template.svg"])

