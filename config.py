import configparser
import yaml
import os, sys, re
import copy
import logging

class ConfigClass:
    def __init__(self):
        self.init_done = False
        self.cwd = os.getcwd() + '/'
        self.config = {}

        self.parser = configparser.ConfigParser(comment_prefixes=("#","!")) 

        yaml_file = open(os.path.join(os.path.dirname(__file__), 'config.yaml'))
        self.config_defaults = yaml.load(yaml_file, Loader=yaml.BaseLoader)['options']
        yaml_file.close()
        self.config_defaults = {k.lower(): v for k, v in self.config_defaults.items()}

        # Pass default settings to parser
        for key in self.config_defaults:
            kattr = self.config_defaults[key]
            self.parser['DEFAULT'][key.lower()] = kattr['default']

    # Check legality and set values
    def setValues(self, data_type, vname, values=None):
        config_error = False
        if data_type == 'integerlist':
            try:
                self.config[vname] = self.parser.get("DEFAULT",vname)
                self.config[vname] = tuple([int(field) for field in self.config[vname].split()])
            except:
                config_error = True
                sys.stderr.write('Option "%s" should be a list of integers\n' % (vname))
                pass
        if data_type == 'floatlist':
            try:
                self.config[vname] = self.parser.get("DEFAULT",vname)
                self.config[vname] = tuple([float(field) for field in self.config[vname].split()])
            except:
                config_error = True
                print(self.config[vname].split())
                print([float(field) for field in self.config[vname].split()])
                sys.stderr.write('Option "%s" should be a list of floats\n' % (vname))
                pass
        elif data_type == 'integer':
            try:
                self.config[vname] = self.parser.getint('DEFAULT',vname)
            
            except:
                config_error = True
                sys.stderr.write('Option "%s" should be a integer\n' % (vname))
        elif data_type == 'float':
            try:
                self.config[vname] = self.parser.getfloat('DEFAULT', vname)
            except:
                config_error = True
                sys.stderr.write('Option "%s" should be a float\n' % (vname))
        elif data_type == 'boolean':
            try:
                self.config[vname] = self.parser.getboolean('DEFAULT',vname)
            except:
                config_error = True
                sys.stderr.write('Option "%s" should be a boolean\n' % (vname))
                pass
        elif data_type == 'string':
            if values is not None:
                pvalue = self.parser.get("DEFAULT",vname)
                if pvalue not in values:
                    config_error = True
                    sys.stderr.write('Option "%s" should be one of: %s\n' % (vname, ", ".join(values)))
                else:
                    self.config[vname] = pvalue
            else:
                self.config[vname] = self.parser.get("DEFAULT",vname)
        return config_error


    def initialize(self, config_file="INCAR"):
        if os.path.isfile(config_file):
            with open(config_file) as stream:
                lines = self.fix_expr(stream.read())
                self.parser.read_string("[DEFAULT]\n"+lines)
            self.config_path = os.path.abspath(config_file)
        else:
            print("Specified INCAR '%s' does not exist" % ''.join(os.path.abspath(config_file)), sys.stderr)
            sys.exit(2)

        
        # Check that all options in config.ini are in the configparser
        config_err = False
        b = self.parser.defaults()

        # What happens if INCAR contains unsupported options?
        section_diff = list(set(b) - set(self.config_defaults.keys()))
        if len(section_diff) > 0:
            config_err = True
            raise NameError('unknown option "%s" \n' % (", ".join(section_diff)))

        # Check type of input settings and assignments
        # for section in psections:
        defaults = self.config_defaults
        for op in b:
            if 'values' in defaults[op]:
                config_err = self.setValues(defaults[op]['kind'], op, values=defaults[op]['values'])
            else:
                config_err = self.setValues(defaults[op]['kind'], op)
        if config_err:
            sys.exit(2)

    def fix_expr(self,line):
        # this is the only one I've seen used in vasp.
        # if there are others, we can certainly implement
        # not sure how to handle the multiline strings . . . (right now)
        lines = line.replace("; ","\n").replace(";","\n") # VASP allows ; to write on same line
        lines = lines.replace(" \\\n"," ") # VASP allows multiline arrays
        lines = re.sub(r"[ \t]{2,}"," ",lines) # fix spacing that resulted from previous linef
        # replace NUM*NUM with what it evaluates to
        iter1 = re.finditer(r'\s[0-9]+\*[0-9]+\s',lines)
        for m in iter1:
            begin,end = m.span()
            idx = lines[begin:end].find("*")+begin
            sub = f" {' '.join(int(lines[begin:idx])*[lines[idx+1:end]])} "
            lines = lines[:begin] + sub + lines[end+1:]
        return lines

if __name__ == "__main__":
    config = ConfigClass()
    config.initialize()
    print("Parameters are stored in a dictionary:")
    print("  ",config.config)
    print("    To fetch the particular parameter, use 'config.config[<nameOFparameter>]'")

