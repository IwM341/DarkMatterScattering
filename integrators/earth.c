inline double earth_profile(const double r)
// returns Earth's density profile in ton/m^3=g/cm^3
// r - in km
// from 1205.0930
// SD 08.07.2012
{
  const double rr = r/6371.0;
  double dens;
  if(r < 1217.1)
    {
      dens = 13.0122 - 8.45292*rr*rr;
    }
  else
    {
      if(r < 3458.7)
	{
	  dens = 12.5842 - 1.69929*rr - 1.94128*rr*rr
	    - 7.11215*rr*rr*rr;
	}
      else
	{
	  if(r < 5701.0)
	    {
	      dens = 6.8143 - 1.66273*rr - 1.18531*rr*rr;
	    }
	  else
	    {
	      if(r < 5951.0)
		{
		  dens = 11.1198 - 7.87054*rr;
		}
	      else
		{
		  if(r < 6336.0)
		    {
		      dens = 7.15855 - 3.85999*rr;
		    }
		  else
		    {
		      if(r < 6351.0)
			{
			  dens = 2.92;
			}
		      else // rr < 6371.0
			{
			  dens = 2.72;
			}
		    }
		}
	    }
	}
    }
  return dens; // in g/cm^3 = ton/m^3

}
