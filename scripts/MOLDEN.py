

class Read_molden(RawFile):
    def __init__(self, filename):
        super(Read_molden, self).__init__(filename)
        self.pt = PeriodicTable()
    
    def findfreq(self):
        startline = self.searchForString('[FREQ]', False)
        if startline == -1:
            raise KeyError('No FREQ attribute on molden file!')
        
        num = int(self.fileContents[startline].split()[1]) - 6
        self.freq = []
        for i in range(num):
            self.freq.append( float( self.fileContents[startline+1+i].split()[0] ) )
        
        self.freq = np.array(self.freq)
        return self.freq
    
    def findNormalMode(self):
        startline = self.searchForString('[FR-NORM-COORD]', False)
        if startline == -1:
            raise KeyError('No FR-NORM-COORD attribute on molden file!')
        
        atom_num = int(self.fileContents[startline].split()[1])//3
        mode_num = atom_num*3 - 6
        self.normalMode = []
        for i in range(mode_num):
            vib = []
            self.normalMode.append(vib)
            for j in range(atom_num):
                vib.extend(self.fileContents[startline+i*(atom_num+1)+j+2].split())
        
        self.normalMode = np.array(self.normalMode, dtype=float)
        return self.normalMode
    
    def findoptGeom(self, frame=-1):
        '''
        frame specify which frame as optimized geometry
        default -1 (the final one)
        '''
        startline = self.searchForString('[GEOMETRIES]', False)
        if startline == -1:
            raise KeyError('No GEOMETRIES attribute on molden file!')
        
        num = int(self.fileContents[startline+1].split()[0])
        # frameline = startline+1
        if frame == -1:
            self.opt_xyz = self.fileContents[-num-2:]
            # return self.opt_xyz
            return ''.join(self.opt_xyz)
    
    def get_mass(self):
       return self.mass
        
    
    def get_opt_cord(self):
        try: self.opt_xyz
        except: self.findoptGeom()

        self.elements = []
        self.mass = []
        self.opt_position = []

        for line in self.opt_xyz[2:]:
            i = line.split()
            try:
                Z = self.pt.getZ(i[0])
            except:
                Z=int(i[0])
            self.elements.append(Z)
            self.mass.append([self.pt.mass[Z]]*3)
            self.opt_position.extend(i[1:4])
        
        self.mass = np.array(self.mass).reshape(-1)
        # n = int(self.opt_xyz[0].split()[0])
        # for i in range(len(self.opt_xyz)):
        #     if i%(n+2)==0:
        #         p = []
        #     elif i%(n+2)==1: continue
        #     else:
        #         p.extend(self.opt_xyz[i].split()[1:4])
        return np.array(self.opt_position, dtype=np.float64)        
    
    def optGeom2xyz(self, xyzName):
        with open(xyzName, 'w') as f:
            try: self.opt_xyz
            except: self.findoptGeom()
            f.writelines(self.opt_xyz)