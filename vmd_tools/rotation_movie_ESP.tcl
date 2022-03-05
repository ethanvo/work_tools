proc make_rotation_movie_files {} {
	set frame 0
	for {set i 0} {$i < 360} {incr i 1} {
		set filename frames/espmovie.[format "%04d" $frame].tga
		render TachyonInternal $filename
		incr frame
		rotate y by 1
	}
}
