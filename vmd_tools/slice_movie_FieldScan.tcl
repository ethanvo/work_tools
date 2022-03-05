proc make_slice_movie_files {} {
	set frame 0
	for {set i 0} {$i < 250} {incr i 1} {
		set filename frames/fieldscanmovie.[format "%04d" $frame].tga
		mol modstyle 2 0 VolumeSlice [expr $i * 0.004] 2.000000 2.000000 2.000000
		render TachyonInternal $filename
		incr frame
	}
}