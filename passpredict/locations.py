from orbit_predictor.locations import Location as LocationBase


class Location(LocationBase):

    @property
    def lat(self):
        return self.latitude_deg

    @property
    def lon(self):
        return self.longitude_deg

    @property
    def h(self):
        return self.elevation_m

    def __repr__(self):
        s = '<Location '
        if self.name:
            s += self.name + ' '
        s += f'({self.latitude_deg}, {self.longitude_deg})'
        s += '>'
        return s
