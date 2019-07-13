#-*- coding:utf-8 -*-

class IOFile:
    """
    A class to input & output file
    The file name with .bam/.gz can be recognized
    """
    def __init__(self, fn, mode='r'):
        assert mode in ('r', 'w')
        self.mode = mode
        self.fn = fn
        if mode == 'r':
            if fn.endswith('.gz'):
                import gzip
                self.handler = gzip.open(fn, 'rt')
                self.filetype = 'GZ'
            elif fn.endswith('.bam'):
                self.handler = os.popen("samtools view -h "+fn, 'r')
                self.filetype = 'BAM'
            else:
                self.handler = open(fn, 'rt')
                self.filetype = 'NORMAL'
        elif mode == 'w':
            if fn.endswith('.gz'):
                import gzip
                self.handler = gzip.open(fn, 'wt')
                self.filetype = 'GZ'
            elif fn.endswith('.bam'):
                import tempfile
                self.tmp_file = tempfile.mkstemp(suffix=".tmp.sam")[1]
                self.handler = open(self.tmp_file, 'wt')
                self.filetype = 'BAM'
            else:
                self.handler = open(fn, 'wt')
                self.filetype = 'NORMAL'
        else:
            print("Error")
    
    def writelines(self, context):
        return self.handler.writelines(context)
    
    def readline(self):
        return self.handler.readline()
    
    def readlines(self):
        return self.handler.readlines()
    
    def close(self):
        self.handler.close()
        if self.mode == 'w' and self.filetype == "BAM":
            os.system("samtools view -bh {} > {}".format(self.tmp_file, self.fn))
            os.system("rm "+self.tmp_file)
    
    def endOfFile(self):
        if self.handler.closed:
            return True
        current = self.handler.tell()
        line = self.handler.readline()
        self.handler.seek(current)
        return True if line=="" else False
    
    def __del__(self):
        if hasattr(self, 'tmp_file') and os.path.exists(self.tmp_file):
            os.system("rm "+self.tmp_file)

