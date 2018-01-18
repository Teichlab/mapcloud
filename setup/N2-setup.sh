#spin up a mapcloud from the snapshot made in the previous step

#log onto delta.internal.sanger.ac.uk, instances, launch instance, name it something informative
#instance flavour: k1.2xlarge (28 cores) - best use of 116 available cores and 4 allocatable floating IPs
#boot from image, set image to the snapshot you created in the previous step
#access & security tab, tick default + ssh + icmp
#networks tab, make sure cloudforms is dragged in
#once spawned, press the little arrow on the far right of the instance's row and associate a floating IP
#volumes, create, name it something informative, size 2048gb
#once spawned, press the little arrow on the far right of the volume's row and attach to the created instance

#setting up the mount in /mnt
sudo mkfs.ext4 /dev/vdb
sudo mount /dev/vdb /mnt
#the problematic mount file from earlier also respawns when a new cloud is made
#so let's swap it to point at the mount we just created
sudo sed 's/data1/mnt/g' -i /etc/fstab

#move ownership of the mounted drive to ubuntu and "jog" it to see it works
#if something is off, the jog will hang badly
sudo chown -R ubuntu: /mnt
cd /mnt && dd if=/dev/zero of=deleteme oflag=direct bs=1M count=1024 && rm deleteme

#acquire the code to do things with!
cd /mnt && git clone https://github.com/Teichlab/mapcloud

#need to configure the internal.sanger.ac.uk thing again, this time more aggressively
#as things go heavily south for some reason when snapshotting happens
echo 'supersede domain-name "internal.sanger.ac.uk";' | sudo tee -a /etc/dhcp/dhclient.conf
sudo dhclient -r
sudo reboot & ( sleep 30; echo 'b' > /proc/sysrq-trigger )
#this restarts the thing. ssh back in
