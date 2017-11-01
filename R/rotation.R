`Rotation` <-
function(xy, angle)
{
	xy<-as.matrix(xy)
	### Find cos and sin of the angle
	cos.angle <- cos(angle)
	sin.angle <- sin(angle)
	
	### Rotate the set of coordinates
	xy.rot<-xy %*% t( matrix(c(cos.angle,sin.angle, -sin.angle, cos.angle), 2,2) )
	
	return(xy.rot)
}

