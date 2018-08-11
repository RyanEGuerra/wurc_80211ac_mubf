# wurc_80211ac_mubf
Library and examples to generate, transmit, and decode IEEE 802.11ac and IEEE 802.11af standard Multi-User MIMO (MU-MIMO) packets of arbitrary dimension and channel bandwidth using the WARPv3 software-defined radio framework from Mango Communications using MATLAB scripts. Optionally, a software switch will allow you to perform the same experiments using a WARPv3 kit and a WURC kit from Skylark Wireless to transmit 802.11af MU-MIMO packets over UHF channels.

This code will form 802.11ac sounding packets, transmit them over the air to a number of receivers to estimate the channel, then generate a multi-user transmission of arbitrary size and channel bandwidth using zero-forcing beamforming to a set of those receivers. The data is decoded and reception metrics are provided as output.

This framework can be used in a real environment to test actual MUBF performance under various conditions. It also has convenient spots to insert your own beamforming user selection or beamforming algorithms, provided to you stick to the 802.11ac/af frame format. I'm releasing this open-source in the hopes that it will assist anyone looking to learn about wireless systems, beamforming, or software-defined radios.

## Known Bugs
1. 4x3 transmissions have a bug in the beamforming calculation logic somewhere. I recommend sticking to 4x2 or 4x4 transmissions.

## How do I start?
I recommend you start at the top level debug_* scripts that provide known working implementations of the library that you can use as a starting point.

You should be able to set the variable USE_OTA_HARDWARE=0 to avoid initializing and transmitting over a physical radio channel with hardware in the loop. Instead, packets are transmitted through a virtual channel, but everything from encoding to decoding is performed exactly the same. This would be a great way to set up simulations or experiments without actual radio hardware or under simulated channel conditions.

This will also make sure that you have you MATLAB environment set up properly. Adding hardware in the loop requires advanced knowledge to configure and set up complex software defined radiom systems. At this time, I am NOT able to assist you with that task.

## WARPLab Versions
Images of the compatible WARPLab image and SD card flashing scripts are available within the /warp_configs directory. This library uses WARPLab version 7.4.0. No other versions have been tested.

## History
I developed this code while working on my Ph.D. thesis. I used it to test hardware systems and wireless protocols that were being developed at that time and to run small tests and experiments in the lab. It was also used by other researchers with some customizations that have not yet been published. Quite frankly, I was often surprised at how well it works and I hope others can learn something from it.

## License
This all is released under the Apache-2.0 license. Best of luck, everyone! Go nuts!
