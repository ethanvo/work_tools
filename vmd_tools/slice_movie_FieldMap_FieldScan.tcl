proc make_slice_movie_files {} {
	set frame 315
	for {set i 0} {$i < 250} {incr i 1} {
		set filename frames/fieldmapscanmovie.[format "%04d" $frame].tga
		mol modstyle 3 0 VolumeSlice [expr $i * 0.004] 2.000000 2.000000 2.000000
		render TachyonInternal $filename
		incr frame
	}
}