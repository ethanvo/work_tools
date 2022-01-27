proc tricolor_scale {} {
  set color_start [colorinfo num]
  display update off
  for {set i 0} {$i < 1024} {incr i} {
    set r [expr {1.0 - ($i * 0.5 / 1023.0)}];  set g [expr {1.0 - ($i / 1023.0)}];  set b [expr {1.0 - ($i * 0.5 / 1023.0)}]
    color change rgb [expr $i + $color_start     ] $r $g $b
  }
  display update on
}
