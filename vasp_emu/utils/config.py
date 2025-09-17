import configparser
import yaml
import os, sys, re
import copy
import logging, warnings

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
        #TODO: If NSW=-1 or 0 and IBRION is unset, the default value should be -1 to indicate no optimizer
        #TODO: ISIF=0 for default if IBRION=0 or LHFCALC = .TRUE.
        # Pass default settings to parser
        # SMASS: Only allows -3, -2, -1, and floats greater than or equal to 0
        # TEEND: while the default says 0.0, it will be set to TEBEG if unspecified
        # ANdersen_prob must be between 0 and 1
        for key in self.config_defaults:
            kattr = self.config_defaults[key]
            self.parser['DEFAULT'][key.lower()] = kattr['default']
    # Check legality and set values
    def setValues(self, data_type, vname, values=None, rules=None):
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
            else:
                if rules is not None:
                    for rule in rules:
                        config_error = self.check_valid(vname,rule) or config_error
        elif data_type == 'float':
            try:
                self.config[vname] = self.parser.getfloat('DEFAULT', vname)
            except:
                config_error = True
                sys.stderr.write('Option "%s" should be a float\n' % (vname))
            else:
                if rules is not None:
                    for rule in rules:
                        config_error = self.check_valid(vname,rule) or config_error
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
        b = self.parser.defaults() # it's not actually defaults, just easiest way to retrieve
        section_diff = list(set(b) - set(self.config_defaults.keys()))
        if len(section_diff) > 0:
            config_err = True
            warnings.warn('unknown option(s): "%s" \n' % (", ".join(section_diff)),stacklevel=2) # is a warning necessary?
            # What happens if INCAR contains unsupported options?
            for val in section_diff: del b[val]

        # Check type of input settings and assignments
        # for section in psections:
        defaults = self.config_defaults
        for op in b:
            if 'values' in defaults[op] and 'rules' in defaults[op]:
                err = self.setValues(defaults[op]['kind'], op, values=defaults[op]['values'], rules=defaults[op]['rules'])
            elif 'values' in defaults[op]:
                err = self.setValues(defaults[op]['kind'], op, values=defaults[op]['values'])
            elif 'rules' in defaults[op]:
                err = self.setValues(defaults[op]['kind'], op, rules=defaults[op]['rules'])
            else:
                err = self.setValues(defaults[op]['kind'], op)
            # this means that we won't accidently reset config_err to F
            # but still allows us to print all the values that are wrong before exiting
            if not config_err and err:
                config_err = True
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

    # Usually, should only be used for ints and floats that have a wide range of allowed values
    def check_valid(self,vname,rule):
        config_err = True
        txt = ""
        #TODO: Maybe  also AND?, depends
        if ".OR." in rule:
            sub_rules = rule.split(".OR.")
        else:
            sub_rules = [rule]

        for i,r in enumerate(sub_rules):
            comp = r.strip()[:4]
            val = r.strip()[4:].strip()
            if comp == ".GE.":
                err = not self.config[vname] >= type(self.config[vname])(val)
            elif comp == ".GT.":
                err = not self.config[vname] > type(self.config[vname])(val)
            elif comp == ".LE.":
                err = not self.config[vname] <= type(self.config[vname])(val)
            elif comp == ".LT.":
                err = not self.config[vname] > type(self.config[vname])(val)
            elif comp == ".EQ.":
                err = not self.config[vname] == type(self.config[vname])(val)
            elif comp == ".NE.":
                err = not self.config[vname] != type(self.config[vname])(val)
            else:
                raise RuntimeError("Unknown rule given:",rule)
            config_err = config_err and err
            txt += f"({vname.upper()} {r.strip()})" + (" .OR. " if i < len(sub_rules)-1 else "")
        if config_err:
            sys.stderr.write(f"Option '{vname.upper()}' does not meet the following criteria: ")
            sys.stderr.write(f"{txt}\n")
        return config_err

if __name__ == "__main__":
    config = ConfigClass()
    config.initialize()
    print("Parameters are stored in a dictionary:")
    print("  ",config.config)
    print("    To fetch the particular parameter, use 'config.config[<nameOFparameter>]'")

