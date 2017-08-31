#ifndef __CPPLibs_MGConfig_H__
#define __CPPLibs_MGConfig_H__


namespace MGConfig {


class JobOpt {
    public :
        JobOpt(int argc, char* argv[]);
        ~JobOpt() {}
        
        const bool& operator() () const { return succ_; }

        const std::string& type() const { return type_; }
        const std::string& ipath() const { return ipath_; }
        const std::string& opath() const { return opath_; }
        
        const long& gi() const { return gi_; }
        const long& gs() const { return gs_; }

        const std::vector<std::string> flist() const { return flist_; }

    private :
        bool        succ_;
        std::string type_;
        std::string ipath_;
        std::string opath_;
        long        gi_;
        long        gs_;
        
        std::vector<std::string> flist_;
};


} // namesapce MGConfig





/*

// 4-> : Title variable type value


namespace MgntConfig {
	typedef std::string   STRING;
	typedef long long int INTEGRAL;
	typedef double        FLOAT;
	union TYPE {
		//std::string   STR;
		long long int INT;
		double        FLT;
	};

	class VariantBase;
	class Config;
}

class MgntConfig::VariantBase {
	public :
		VariantBase() {}
		~VariantBase() {}
		
		inline void setType(const char type) { fType = type; }
		inline const char getType() { return fType; }
		inline MgntConfig::TYPE& operator()() { return fMember; }

	protected :
		char fType;
		MgntConfig::TYPE fMember;
};

// TODO : multi type value
class MgntConfig::Config {
	public :


	public :
		Config() {}
		~Config() {}

		void Register(const std::string& title, const std::vector<std::string>& keys) {}
		void Register(const std::string& title, const std::string& key, const std::string& val = "") {}

		//std::map Get(const std::string& title);
		//std::string& Get(const std::string& title, std::string& key);

	protected :
		//std::map<std::string, std::map<std::string, std::string> > > fConfig;
		std::map<std::string, VariantBase& > fConfig;
};

*/
#endif // __CPPLibs_MGConfig_H__
